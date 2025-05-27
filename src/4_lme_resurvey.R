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
})

# load ReSurveyEurope data
d.initial <- './data/ReSurveyEU_clean.csv.xz' %>%
  read_csv(show_col_types = F)

# set ind names to analyze
ind.names <- c('EIV_M', 'EIV_N', 'EIV_T', 'EIV_L', 'EIV_R')

# set minimum proportion of species with available EIV to include a plot in the analyses
treshold.of.EIVE.species = 0.8

# count number of resurveys involved
d.initial %>%
  select(resurv_id, contains('n.')) %>%
  gather('k', 'v', contains('n.')) %>%
  filter(v >= treshold.of.EIVE.species) %>%
  pull(resurv_id) %>%
  unique() %>%
  length()

# count number of plots involved
d.initial %>%
  select(plot_id, contains('n.')) %>%
  gather('k', 'v', contains('n.')) %>%
  filter(v >= treshold.of.EIVE.species) %>%
  pull(plot_id) %>%
  unique() %>%
  length()

# temporal range of observations 
d.initial$year %>% range

#### 2. Run LME for each EIV variable ####

res <- list()
for (ind.name in ind.names) {
  
  #sti = Sys.time()
  
  # copy initial dataset
  d = d.initial
  
  # Set name of the response variable
  names(d)[which(names(d) %in% paste0('cm.', ind.name))] <- 'eiv'
  names(d)[which(names(d) %in% paste0('n.', ind.name))]  <- 'th'
  
  # apply filters
  d <- d %>%
    # filter based on proportion of species with EIV values available
    filter(th >= treshold.of.EIVE.species) %>%
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
  
  # fit the LME model
  lme_model <- lmer(eiv ~ habitat * time + plot_size_log + (1 |dataset_id / resurv_id),
                    data = d)

  # get slope for time and habitat
  suppressMessages(
    habitat_yearslopes <- emtrends(lme_model, specs = ~ habitat, var = 'time')
  )
  habitat_yearslopes
  
  res[[ind.name]] <- habitat_yearslopes
  
  #print(Sys.time() - sti)
}

#### 3. Plot LME coefficients ####
# aggregate model results
res.dat <-
  res %>% 
  map(as.data.frame) %>% 
  bind_rows(.id = 'ind.name') %>%
  left_join(data.frame(
    ind.name = c('EIV_M', 'EIV_N', 'EIV_T', 'EIV_L', 'EIV_R'),
    eiv_name = c('Moisture','Nutrients','Temperature','Light','Soil reaction'), 
    by = 'ind.name'))

# transform EIV variables to factor
res.dat$eiv_name <-
  factor(res.dat$eiv_name, rev(c('Moisture','Nutrients','Temperature','Light','Soil reaction')))

# visualize results
p <- ggplot(res.dat, aes(
  x=eiv_name, 
  y=time.trend, # EIV change per decade
  col=eiv_name)
  ) + 
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_pointrange(aes(ymin=asymp.LCL, ymax=asymp.UCL), fatten = 3) +
  facet_wrap(~habitat
             #, scales='free'
  ) +
  coord_flip() +
  labs(y=expression(paste(Delta, " CM", ""[EIV], " per decade")), x='') +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    legend.position = 'none',
    axis.text = element_text(color = 'black', size=13),
    axis.title = element_text(size = 13),
    strip.text.x = element_text(hjust = 0, margin=margin(l=0), size=15, face=2),
    strip.background = element_blank()
  ) +
  scale_color_manual(values=RColorBrewer::brewer.pal(5,'Dark2'))

# folder where to store figures 
pth2fig <- './fig/'

# export
ggsave('ReSurv.LME.change.svg', 
       p, 
       path = pth2fig,
       width = 7, height = 7)

# quit
quit(save = 'no')