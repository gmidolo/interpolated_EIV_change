################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Interpolate over time (predict models at each site for each year)

################################################################################

#### 1. Set up data ####

# standardize plot size to median by habitat while predicting EIV values?
standardize_plot.size = TRUE

# folder where model results are stored
pth2export <- './models/'

# folder where data with predictions (in .csv) will be exported
pth2preds <- './preds/'

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(matrixStats) # used for weighted-median plot size calculation
})

# source help functions 
source('./src/0_helpfunctions.R')

# get initial data set
set.seed(123)
dat.initial <-
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope( # Load ReSurveyEU (static, using one random point in the survey)
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
    )[['traintest_data']]  
  )

# calculate habitat plot size medians
if (standardize_plot.size) {
  hab_plot.size <- dat.initial %>%
    group_by(dataset, habitat) %>%
    summarise(plot_size_median_in_db = median(plot_size),
              n_in_db = n()) %>%
    ungroup() %>%
    group_by(habitat) %>%
    summarise(plot_size_median = weightedMedian(x = plot_size_median_in_db, w =
                                                  n_in_db))
}
glimpse(hab_plot.size)

#### 2. Interpolate EIV change metrics ####
# set up response variable names
ind.names <- c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R')

# start interpolating (predicting) EIVs
for(ind.name in ind.names){
  
  # start message
  print(paste0('Start interpolation for ', ind.name))

  # define start time
  st = Sys.time()
  
  # copy initial data set
  dat <- dat.initial
  
  # rename focal variable
  names(dat)[which(names(dat) == paste0('n.', ind.name))] <- 'treshold'
  names(dat)[which(names(dat) == paste0('cm.', ind.name))] <- 'eiv'

  # apply filters accordingly to the raw data
  treshold.of.EIVE.species = 0.8 
  dat <- dat %>%
      filter(treshold >= treshold.of.EIVE.species)

  # final selection of variables needed for modeling
  dat <- dat %>%
      select(plot_id, eiv, x, y, elev, year, habitat, n, plot_size)

  # path to focal model
  pth_model <-
    list.files(pth2export,
               pattern = paste0('RF.last_fit_', ind.name),
               full.names = T)
  
  # load model
  m <- pth_model %>% read_rds()

  # extract workflow for predictions
  m_wf <- extract_workflow(m)

  # define years for prediction
  pred_years <- seq(1960, 2020, 1)

  # copy data for predictions
  pred_dat <- dat 
  
  ## predict over on a subset of plots only?
  # pred_dat <- dat %>% sample_n(1000)
  
  # start predicting for each year
  for (i in pred_years) {
    # cat(i, ' - ')
    pdi <- pred_dat
    pdi <- pdi %>% mutate(year = i) # change year
    if (standardize_plot.size) { # standardize plot size to median by habitat?
      pdi <- pdi %>%
        select(-plot_size) %>% # substitute fixed plot size for a given habitat
        left_join(hab_plot.size, 'habitat') %>%
        rename(plot_size = plot_size_median)
    }
    # get predictions on the new data
    pred_dat[paste0('eiv_pred_', i)] <- predict(object = m_wf, new_data = pdi) %>% pull(.pred)
  }
  
  ## Calculate metrics of change for each plot
  
  # calculate changes for focal time periods: 1960, 1980, 2000, 2020
  pred_focalchange <- pred_dat %>%
    select(plot_id, contains('_pred_')) %>%
  # Calculate percentages of changes:
    mutate(eiv_perc.change_1960.1980 = 100*((eiv_pred_1980 - eiv_pred_1960)/eiv_pred_1960)) %>%
    mutate(eiv_perc.change_1980.2000 = 100*((eiv_pred_2000 - eiv_pred_1980)/eiv_pred_1980)) %>%
    mutate(eiv_perc.change_2000.2020 = 100*((eiv_pred_2020 - eiv_pred_2000)/eiv_pred_2000)) %>%
    mutate(eiv_perc.change_1960.2020 = 100*((eiv_pred_2020 - eiv_pred_1960)/eiv_pred_1960)) %>%
  # Calculate log response ratios:      
    mutate(eiv_lnRR.change_1960.1980 = log(eiv_pred_1980/eiv_pred_1960)) %>%
    mutate(eiv_lnRR.change_1980.2000 = log(eiv_pred_2000/eiv_pred_1980)) %>%
    mutate(eiv_lnRR.change_2000.2020 = log(eiv_pred_2020/eiv_pred_2000)) %>%
    mutate(eiv_lnRR.change_1960.2020 = log(eiv_pred_2020/eiv_pred_1960)) %>%
  # Calculate absolute change:      
    mutate(eiv_abs.change_1960.1980 = eiv_pred_1980-eiv_pred_1960) %>%
    mutate(eiv_abs.change_1980.2000 = eiv_pred_2000-eiv_pred_1980) %>%
    mutate(eiv_abs.change_2000.2020 = eiv_pred_2020-eiv_pred_2000) %>%
    mutate(eiv_abs.change_1960.2020 = eiv_pred_2020-eiv_pred_1960) %>%
  # round to 3 digits
    mutate_all(round, 3) %>%
    # remove preds
    select(-contains('_pred_'))

  # calculate linear slope estimates, estimated with predicted S ~ year (period 1960-2020)
  pred_lmslope <- pred_dat %>%
    select(plot_id, contains('_pred_')) %>%
    gather('year', 'eiv_pred', contains('_pred_')) %>%
    mutate(year = as.numeric(gsub("\\D", "", year))) %>%
    group_by(plot_id) %>%
    do(tidy(lm(eiv_pred ~ year, data = .))) %>%
    filter(term == 'year') %>%
    select(plot_id, estimate, std.error) %>%
    setNames(c('plot_id', 'lm.slope_estimate', 'lm.slope_std.error'))

  # Calculate linear slope estimates, estimated with predicted S ~ year but for each 20 yrs intervals
  period = c('1960.1980','1980.2000','2000.2020')
  min_p = c(1960, 1980, 2000)
  max_p = c(1980, 2000, 2020)
  res_lm_periods = list()
  for(i in period) {
    # cat('PERIOD: ', i, ' - ')
    res_lm_periods[[i]] <- pred_dat %>%
      select(plot_id, contains('_pred_')) %>%
      gather('year', 'eiv_pred', contains('_pred_')) %>%
      mutate(year = as.numeric(gsub("\\D", "", year))) %>%
      filter(year >= min_p[which(period %in% i)]) %>%
      filter(year <= max_p[which(period %in% i)]) %>%
      group_by(plot_id) %>%
      do(tidy(lm(eiv_pred ~ year, data = .))) %>%
      filter(term == 'year') %>%
      select(plot_id, estimate, std.error) %>%
      setNames(c('plot_id', paste0('lm.slope_estimate',i), paste0('lm.slope_std.error',i)))
  }

  # Join predictions made
  pred_dat <- pred_dat %>%
              left_join(pred_focalchange, 'plot_id') %>%
              left_join(pred_lmslope, 'plot_id') %>%
              left_join(res_lm_periods$`1960.1980`, 'plot_id') %>%
              left_join(res_lm_periods$`1980.2000`, 'plot_id') %>%
              left_join(res_lm_periods$`2000.2020`, 'plot_id')

  # add EIV name column
  pred_dat <- pred_dat %>% mutate(eiv_name=ind.name, .after=plot_id)

  # export predictions in .csv
  pred_dat %>%
    write_csv(paste0(pth2preds, ind.name, '.preds.rf.csv.gz'))
  
  # measure time duration
  ett <- Sys.time()-st
  print(paste('Interpolation for ', ind.name, ' done: ', round(ett[[1]], 2),  units(ett), 'to complete'))

}

# quit
quit(save='no')

# # final summary (optional): count 
# ind.names <- c('EIV_M','EIV_N','EIV_T','EIV_L','EIV_R')
# plots_in_intrpl <- ind.names %>%
#   map(\(x){
#     read_csv(paste0(pth2preds, ind.name, '.preds.rf.csv.gz')) %>% 
#       select(plot_id,year) %>% 
#       unique()
#   }) %>%
#   bind_rows() %>%
#   unique()
# 
# plots_in_intrpl$year %>% range
# nrow(plots_in_intrpl[plots_in_intrpl$year <= 2020 & plots_in_intrpl$year>=1960,]) 