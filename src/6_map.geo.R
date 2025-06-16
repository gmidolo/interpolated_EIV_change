################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Plot geographic maps of EIV changes across Europe

################################################################################


#### 1. Prepare data ####

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(rnaturalearth)
})

# set indicator names
ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

# source help functions 
source('./src/0_helpfunctions.R')

# folder where data with predictions (in .csv) were exported
pth2preds <- './preds/'

# folder where to store figures and diagnostic output
pth2fig <- './fig/'

# prepare data for interpolation
set.seed(123)
dat.initial <-
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope( # Load ReSurveyEU (static, using one random point in the survey)
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
    )[['traintest_data']]  
  )

# load predictions for each EIV variable
dat <- ind.names$eiv_name_raw %>%
  map(\(ind.name){
  # get plot id to discard
  th <- dat.initial[,c('plot_id', paste0('n.', ind.name))] %>%
        setNames(c('plot_id','th')) %>%
        filter(th<0.8)
  # get predictions
  y <- paste0(pth2preds, ind.name, '.preds.rf.csv.gz') %>%
    read_csv(show_col_types = F) %>%
    filter(year >= 1960 & year <= 2020) %>%
    rename(eiv_name_raw = eiv_name) %>%
    left_join(ind.names, 'eiv_name_raw') %>%
    select(-eiv_name_raw) %>%
    select(eiv_name, everything()) %>%
    arrange(plot_id) %>%
    anti_join(th, 'plot_id')
  return(y)
  })

# set new names
names(dat) <- ind.names$eiv_name

# subset needed data
dat <- dat %>%
 map(\(x){
  x %>%
   select(eiv_name, habitat, plot_id, x, y, eiv_abs.change_1960.2020)
 }) %>%
 bind_rows() 

# get EU-countires sf shapes
regions_name <- c('Albania', 'Austria', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria',
                  'Corsica', 'Crete', 'Croatia', 'Czechia', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany',
                  'Greece', 'Hungary', 'Ireland', 'Italy', 'Kosovo', 'Latvia', 'Liechtenstein',
                  'Lithuania', 'Luxembourg', 'Malta', 'Moldova', 'Montenegro', 'Netherlands', 'North Macedonia',
                  'Norway', 'Poland', 'Portugal', 'Romania', 'Sardinia', 'Serbia', 'Sicily',
                  'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine', 'United Kingdom')
bbox_coords <- c(xmin = -1123055, ymin = 3923814, xmax = 2796649, ymax = 8007282)
EU <- ne_countries(scale = 'large', returnclass = 'sf') %>%
  filter(name %in% regions_name) %>%
  st_transform(crs = 25832) %>%
  st_crop(bbox_coords) %>% 
  select(geometry)

# set theme for plotting
my_theme <- theme(
  #panel.background = element_rect(fill = 'white',color='grey70'),
  #panel.grid = element_line(color = 'grey90'),
  #strip.background = element_rect(fill='white', color='black'),
  strip.text = element_text(colour ='black'),
  axis.title=element_blank()
)

#### 2. Rasterize average change per habitat ####

# set base raster
res_km <- 50
r <- rast(res = res_km * 1000,
          extent = ext(EU),
          crs = crs(EU))

# variables to work with, assign time range
var2collect <- data.frame(var = 'eiv_abs.change_1960.2020',
                          min_yr = 1960,
                          max_yr = 2020)

# Funtions used to aggregate metrics of change
mean_at_least_five_plots <- \(x, th = 5) {
  ifelse(length(x) >= th, mean(x), NA)
} # average if there are at least 5 plots
count_at_least_five_plots <- \(x, th = 5) {
  ifelse(length(x) >= th, length(x), NA)
} # average if there are at least 5 plots

dat_i_h <- list()
conta_tot = list()
conta_rast = list()
d2rast <- list()
d2plot <- list()
for (e in ind.names$eiv_name) {
  for (i in var2collect$var) {
    for (h in c('Forest', 'Grassland', 'Scrub', 'Wetland')) {
      
      d_i_h <- dat %>%
        filter(eiv_name == e) %>%
        filter(habitat == h)
      
      dat_i_h[[e]][[h]] <- d_i_h
      
      conta_tot[[e]][[h]] <-
        data.frame(habitat = h, n = nrow(d_i_h))
      
      d2rast[[e]][[h]] <-
        rasterize(d_i_h %>% select(x, y) %>% as.matrix(),
                  r,
                  values = d_i_h[i],
                  fun = mean_at_least_five_plots)
      
      conta_rast[[e]][[h]] <-
        rasterize(
          d_i_h %>% select(x, y) %>% as.matrix(),
          r,
          values = 1:nrow(d_i_h),
          fun = count_at_least_five_plots
        )
      
      d2plot[[e]][[h]] <- d2rast[[e]][[h]] %>%
        as.data.frame(xy = T, na.rm = F) %>%
        rownames_to_column('id') %>%
        rename(mean = values) %>%
        drop_na()
    }
  }
}

d2plot <- d2plot %>%
  map(\(x) {
    bind_rows(x, .id = 'Habitat')
  }) %>%
  bind_rows(.id = 'EIV') 


#### 3. Plotting ####

br = c(-Inf,-1, -0.5, -0.1, 0.1, 0.5, 1, Inf)
lb = c('< 1','-1 - -0.5', '-0.5 - -0.1', '-0.1 - 0.1', '0.1 - 0.5', '0.5 - 1', '> 1')

# prepare data for plotting
d2plot_refined <- d2plot %>%
  # transform continuous var to categories 
  mutate(mean_cat = cut(mean, breaks = br, labels = lb))

# get color palette for each categorical level
cols <-
  hcl.colors(length(levels(d2plot_refined[1, 'mean_cat'])), 'RdYlBu', rev = T)

# simplify EU shape
newEU <- EU %>% st_buffer(1000) %>% st_simplify(dTolerance = 4000)

# convert EIV names to factor
d2plot_refined$EIV <- factor(d2plot_refined$EIV, ind.names$eiv_name)

# plot
p <- ggplot() +
  geom_sf(data=newEU, fill='white', color=NA) +
  geom_raster(data=d2plot_refined, aes(x,y,fill=mean_cat)) +
  facet_grid(EIV ~ Habitat) +
  scale_fill_manual(values = cols) +
  geom_sf(data=newEU, fill=NA, color='grey40', linewidth=.1) +
  labs(fill='Community mean EIV change\n(2020 vs. 1960)')+
  theme(strip.text = element_text(colour ='black', face = 2, size=10),
        axis.title=element_blank(),
        legend.title = element_text(face='bold', size=10),
        legend.position = 'bottom'
  ) +
  guides(fill = guide_legend(reverse=TRUE))

# export
ggsave(paste0(pth2fig, 'EIVchangemap.svg'), p, width = 6, height = 8)

# quit
quit(save = 'no')