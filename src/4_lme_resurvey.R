################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Estimate EIVs trends on time series (ReSurveyEurope) data
# We use linear mixed effect models

################################################################################

#### 1. Prepare data for modelling ####

# Load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(lme4)
  library(emmeans)
  library(rnaturalearth)
})

# load ReSurveyEurope data
d.initial <- './data/ReSurveyEU_clean.csv.xz' %>%
  read_csv(show_col_types = F)

# undesired plots (observations outside the desired range)
outside.period <- d.initial %>%
  filter(year < 1960 | year > 2020)

# apply filters and retain plots that have still more than one survey observation
d.initial <- d.initial %>%
  anti_join(outside.period, by = 'plot_id') %>%
  group_by(resurv_id) %>%
  filter(n() > 1) %>%
  ungroup()

# identify plots (entire time series) with changes in habitat type (e.g. grasslands shifting to scrub)
hab.change.serie <- d.initial %>%
  select(resurv_id, habitat) %>%
  unique() %>%
  group_by(resurv_id) %>%
  summarise(hab_change = n() > 1)
table(hab.change.serie$hab_change)

# apply further filters
d.initial <- d.initial %>%
  # exclude plots (entire time series) with habitat change
  anti_join(hab.change.serie %>% filter(hab_change), 'resurv_id') %>%
  group_by(resurv_id) %>%
  # filter resurvey plot with at least two plots
  filter(n_distinct(plot_id) >= 2) %>%
  # filter resurvey that did not change habitat (EUNIS-ESy lev. 1) between any of the resurvey
  filter(n_distinct(habitat) == 1) %>%
  # make sure plots form the same resurvey id do not belong to different datasets?
  filter(n_distinct(dataset) == 1) %>%
  ungroup()

# set ind names to analyze
ind.names <- c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R')

# set minimum proportion of species with available EIV to include a plot in the analyses
threshold.of.EIVE.species = 0.8

# count number of resurveys involved
d.initial %>%
  select(resurv_id, contains('n.')) %>%
  gather('k', 'v', contains('n.')) %>%
  filter(v >= threshold.of.EIVE.species) %>%
  pull(resurv_id) %>%
  unique() %>%
  length()

# count number of plots involved
d.initial %>%
  select(plot_id, contains('n.')) %>%
  gather('k', 'v', contains('n.')) %>%
  filter(v >= threshold.of.EIVE.species) %>%
  pull(plot_id) %>%
  unique() %>%
  length()

# temporal range of observations 
d.initial$year %>% range

#### 2. Run LME for each EIV variable ####

res <- list()
for (ind.name in ind.names) {
  
  # copy initial dataset
  d = d.initial
  
  # Set name of the response variable
  names(d)[which(names(d) %in% paste0('cm.', ind.name))] <- 'eiv'
  names(d)[which(names(d) %in% paste0('n.', ind.name))]  <- 'th'
  
  # apply filters
  d <- d %>%
    # filter based on proportion of species with EIV values available
    filter(th >= threshold.of.EIVE.species) %>%
    # group by resurvey plot
    group_by(resurv_id) %>%
    # filter resurvey plot with at least two plots
    filter(n_distinct(plot_id) >= 2) %>%
    # filter resurvey that did not change habitat (EUNIS-ESy lev. 1) between any of the resurvey
    filter(n_distinct(habitat) == 1) %>%
    # make sure plots form the same resurvey id do not belong to different datasets?
    filter(n_distinct(dataset) == 1) %>%
    ungroup()  %>%
    # create group indices by dataset name
    group_by(dataset) %>%
    mutate(dataset_id = cur_group_id()) %>%
    ungroup()
  
  # prepare predictors
  d$time <- d$year / 10 # time expressed in decades
  d$plot_size_log <- log10(d$plot_size) # log10 transformed plot size
  d$resurv_id <- as.factor(d$resurv_id) # force resurvey id to factor
  d$dataset_id <- as.factor(d$dataset_id) # force dataset id to factor
  
  # select used data
  d <- d %>%
    select(dataset_id, resurv_id, eiv, habitat, time, plot_size_log)
  
  # number of resurvey plots per habitat
  d.plot.count <- d %>%
    select(resurv_id,habitat) %>%
    unique () %>%
    group_by(habitat) %>%
    summarise(no.plots = n())
  
  # fit the LME model
  lme_model <- lmer(eiv ~ habitat * time + plot_size_log + (1 |dataset_id / resurv_id),
                    data = d)

  # get slope for time and habitat
  suppressMessages(
    habitat_yearslopes <- emtrends(lme_model, specs = ~ habitat, var = 'time')
  )
  habitat_yearslopes
  
  res[[ind.name]] <- habitat_yearslopes
  
}

#### 3. Plot LME coefficients ####
# aggregate model results

res.dat <-
  res %>% 
  map(as.data.frame) %>% 
  bind_rows(.id = 'ind.name') %>%
  left_join(data.frame(
    ind.name = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
    eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction')), 
    by = 'ind.name')

# transform EIV variables to factor
res.dat$eiv_name <-
  factor(res.dat$eiv_name, rev(c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction')))

# visualize results
p <- ggplot(res.dat, aes(
  x=eiv_name, 
  y=time.trend, # EIV change per decade
  col=eiv_name)
  ) + 
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_pointrange(aes(ymin=asymp.LCL, ymax=asymp.UCL), fatten = 4) +
  facet_wrap(~habitat) +
  coord_flip() +
  labs(y = expression(paste(Delta, ' CM', ''[EIV], ' per decade')), x='') +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    legend.position = 'none',
    axis.text = element_text(color = 'black', size=13),
    axis.title = element_text(size = 13),
    strip.text.x = element_text(hjust = 0, margin=margin(l=0), size=15, face=2),
    strip.background = element_blank()
  ) +
  scale_color_manual(values=RColorBrewer::brewer.pal(5,'Dark2')) +
  scale_y_continuous(breaks = c(-0.04, -0.02, 0, 0.02, 0.04, 0.06))
print(p)

# folder where to store figures 
pth2fig <- './fig/'

# export figure
ggsave('ReSurv.LME.change.svg', 
       p, 
       path = pth2fig,
       width = 8, height = 6)

## Plot re-survey location (mini-maps) ##

# get EU-countires sf shapes
EU <- read_rds('./data/EU_shape_map.rds') %>% 
  st_buffer(2000) %>% 
  st_simplify(dTolerance = 4000) %>%
  st_union()

grd_size_km = 75 # grid size
for (i in c('Forest', 'Grassland', 'Scrub', 'Wetland')) {
  # plot occupancy for each habitat
  dip <- d.initial %>%
    group_by(resurv_id) %>%
    filter(n_distinct(plot_id) >= 2) %>%
    filter(n_distinct(habitat) == 1) %>%
    filter(n_distinct(dataset) == 1) %>%
    ungroup() %>%
    filter(habitat == i) %>%
    select(resurv_id, x_mean, y_mean) %>%
    unique()
  dipspat <- st_as_sf(dip,
                      coords = c('x_mean', 'y_mean'),
                      crs = st_crs(EU))
  grd <- st_make_grid(EU, cellsize = c(grd_size_km * 1000, grd_size_km *
                                         1000)) %>%
    st_as_sf()
  dip_grid <- st_filter(grd, dipspat)
  phab <- ggplot() +
    geom_sf(data = EU,
            col = NA,
            fill = '#2d3068') +
    geom_sf(data = dip_grid,
            fill = 'red',
            col = NA) +
    theme_void() +
    theme(legend.position = 'none') +
    theme(panel.background = element_rect(colour = 'black', fill = 'grey85')) +
    ggtitle(paste0('n = ', nrow(dip) %>% prettyNum(big.mark = ",")))
  
  ggsave(
    paste0('ReSurv.minimap.grid.', i, '.pdf'),
    phab,
    path = pth2fig,
    width = 2,
    height = 2
  )
}

# quit
quit(save = 'no')