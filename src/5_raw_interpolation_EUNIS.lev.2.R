################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 17.12.2025
################################################################################

# Description: Summarise interpolation stats using raw (all trees) predictions 
#              across EUNIS-ESy 2 habitat types

################################################################################

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(matrixStats)
})

# standardize plot size to median by habitat?
standardize_plot.size = TRUE

# set indicator names
ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

# source function to format Resurvey data
source('./src/0_helpfunctions.R')

# folder where data with predictions (in .csv) were exported
pth2preds <- './preds/pdp/'

# folder where models trained on EUNIS l. 2 habitats are stored
pth2models <- './models/ESy2models/'

# prepare data for modeling and split train and test dataset
set.seed(123)
dat.initial <-
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope( # Load ReSurveyEU (static, using one random point in the survey)
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
    )[['traintest_data']]  
  )

# standardize plot size to median by habitat while predicting EIV values?
standardize_plot.size = TRUE

# get list of possible ESy2 outcomes
ESy2list <- dat.initial %>% 
  select(ESy) %>% 
  unique %>% 
  mutate(ESy2 = unname(sapply(.$ESy, extract_habitat_lev2))) %>%
  drop_na() %>%
  #remove coastal habitats and palsa mires
  filter(ESy2 != 'N1') %>%
  filter(ESy2 != 'Q3') 

# calc year range for interpolation for each ESy2 habitat
dat.year.ESy2 <- dat.initial %>%
  semi_join(ESy2list, 'ESy') %>%
  left_join(ESy2list, 'ESy')
treshold.of.EIVE.species = 0.8
dat.year.ESy2 <- dat.year.ESy2 %>%
  filter(n.EIV_M >= treshold.of.EIVE.species) %>%
  filter(n.EIV_N >= treshold.of.EIVE.species) %>%
  filter(n.EIV_T >= treshold.of.EIVE.species) %>%
  filter(n.EIV_L >= treshold.of.EIVE.species) %>%
  filter(n.EIV_R >= treshold.of.EIVE.species) 
dat.year.ESy2 <- dat.year.ESy2 %>%
  group_by(ESy2) %>%
  #calculate year range
  summarise(max_yr = round(quantile(year, .95)),
            min_yr = round(quantile(year, .05))) 
head(dat.year.ESy2)

# filter plots based on the interpolation range
dat.year.ESy2.full = list()
for (i in dat.year.ESy2$ESy2) {
  dat.year.ESy2.full[[i]] <- data.frame(
    year = seq(dat.year.ESy2[which(dat.year.ESy2$ESy2 == i), ]$min_yr,
               dat.year.ESy2[which(dat.year.ESy2$ESy2 == i), ]$max_yr,
               1)
  )
}
dat.year.ESy2.full <- dat.year.ESy2.full %>%
  bind_rows(.id = 'ESy2')

# filter out plots sampled outside desired range per ESy2 habitat type
dat.initial <- dat.initial %>%
  left_join(ESy2list, 'ESy') %>%
  semi_join(dat.year.ESy2.full, c('ESy2', 'year'))

# no. of plots
nrow(dat.initial)

#### 2. Predict ####
# obtain average percentage changes per EUNIS-ESy habitat lev.2 and 95%CI
mst = Sys.time()
res <- list() # results of summary stats
for (ind.name in ind.names$eiv_name_raw) {
  st = Sys.time()
  print(paste0('Start interpolation for ', ind.name))
  
  dat <- dat.initial
  
  # assign names to focal EIV variables
  names(dat)[which(names(dat) == paste0('n.', ind.name))] <- 'treshold'
  names(dat)[which(names(dat) == paste0('cm.', ind.name))] <- 'eiv'
  
  # apply filters accordingly to the raw data
  treshold.of.EIVE.species = 0.8
  dat <- dat %>%
    filter(treshold >= treshold.of.EIVE.species)
  
  # final variable selection
  dat <- dat %>%
    select(plot_id, dataset, eiv, x, y, elev, year, ESy2, n, plot_size)
  
  # calculate habitat plot size (weighted) medians
  if (standardize_plot.size) {
    hab_plot.size <- dat %>%
      group_by(dataset, ESy2) %>%
      summarise(plot_size_median_in_db = median(plot_size),
                n_in_db = n()) %>%
      ungroup() %>%
      group_by(ESy2) %>%
      summarise(plot_size_median = weightedMedian(x = plot_size_median_in_db, 
                                                  w = n_in_db))
  }
  
  # load model fit (ESy2-only models)
  pth_model <- paste0(pth2models, 'RF.last_fit_', ind.name, '.rds')
  m <- pth_model %>% read_rds()
  
  # extract workflow for predictions
  m_wf <- extract_workflow(m)
  
  # predict EIVs per year
  year_span <- seq(min(dat$year), max(dat$year), 1)
  pred_res_yr = list()
  for (y in year_span) {
    message(paste0(y, ' - '))
    
    # subset data for prediction in a given year
    pdy <- dat %>%
      # select only ESy2 plots for a given year
      semi_join(dat.year.ESy2.full %>% 
                  filter(year == y) %>% 
                  select(ESy2), by = 'ESy2') %>%
      select(-year)
    pdy$year <- y # replace year for interpolation
    
    # standardize plot size to median by habitat?
    if (standardize_plot.size) {
      pdy <- pdy %>%
        select(-plot_size) %>%
        left_join(hab_plot.size, 'ESy2') %>%
        rename(plot_size = plot_size_median)
    }
    
    # predict
    praw <- predict( m_wf, new_data = pdy, type = 'raw', opts = list(predict.all = T))
    praw <- as_tibble(praw)
    
    # get summary stats
    pred_res_yr[[y]] <- cbind(pdy, praw)  %>%
      group_by(ESy2) %>%
      summarise(across(starts_with('prediction.'), mean)) %>%
      gather('.pred_key', '.pred', starts_with('prediction.')) %>%
      group_by(ESy2) %>%
      summarise(
        median = median(.pred),
        quantile0.05 = quantile(.pred, .05),
        quantile0.95 = quantile(.pred, .95),
        quantile0.025 = quantile(.pred, .025),
        quantile0.975 = quantile(.pred, .975),
        n = n(),
        sd = sd(.pred),
        mean = mean(.pred),
        lower_ci_95 = mean - (1.96 * (sd / sqrt(n))),
        upper_ci_95 = mean + (1.96 * (sd / sqrt(n))),
        .groups = 'drop'
      ) %>%
      mutate(year = y, .before = 1) %>%
      mutate(ind.name = ind.name, .before = 1) %>%
      left_join(pdy %>% group_by(ESy2) %>% summarise(no.plots = n()), 'ESy2') %>%
      left_join(hab_plot.size, 'ESy2')
  }
  
  res[[ind.name]] <- bind_rows(pred_res_yr)
  
  print(Sys.time() - st)
}
print(Sys.time() - mst)

# merge results into a single data frame and export
bind_rows(res) %>% 
  arrange(ind.name, ESy2, year) %>%
  write_csv(paste0(pth2preds, 'ESy2_alltrees.csv'))

# quit
quit(save = 'no')