################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Model validation using time series (ReSurveyEurope) data

################################################################################

#### 1. Data preparation / import ####

# set seed (optional)
set.seed(123) 

# load packages
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(tidymodels)
  }
)

# folder with model results
pth2export <- './models/'

# folder where to store validation results
pth2valid <- './validation/'

# source function to format Resurvey data
source('./src/0_helpfunctions.R')

# Set up response variable names
ind.names <- c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R')

# Set number of iterations
reps = 100

for(ind.name in ind.names){
  # Load tuning results
  tune_res <- list.files(pth2export, pattern = paste0('RF.tune_res_',ind.name), full.names = T) %>%
    read_rds()
  tune_best <-  show_best(tune_res, metric='rmse')[1,]
  # Specify hyperparameters
  hyperpar <- c(trees = 500, min_n = tune_best$min_n, mtry = tune_best$mtry)
  # specify names in data for response variable and threshold proportion of species 
  res.var.name = paste0('cm.', ind.name)
  thr.var.name = paste0('n.', ind.name)
  
  st = Sys.time()
  res = list()
  for(i in 1:reps){ # run the validation for each iteration
    res[[i]] <- run_validation_RF(datatype = 'ReSurv', 
                                  response_var_name = res.var.name,
                                  threshold_value_name = thr.var.name,
                                  threshold_n_min = 0.8,
                                  predictor_formula = '~ x + y + elev + year + habitat + plot_size + n',
                                  training_strategy_ReSurv = 'random',
                                  format_lnRRchange_preds_method = 'end.vs.start',
                                  lnRRchange = FALSE, # obs and pred changes are calculated in absolute EIV change
                                  trees_rf = hyperpar[['trees']], mtry_rf = hyperpar[['mtry']], min_n_rf = hyperpar[['min_n']])
  }
  Sys.time()-st
  
  final_res <- list(
    eval  = res %>% map(\(x){x$eval }) %>% bind_rows(.id = 'rep') 
  )
  
  # export evaluation results only 
  write_csv(final_res$eval, paste0(pth2valid, ind.name,'.ReSu.only.rep_test_validation.csv'))
  
}

# get the stats of the validation results
res <- list()
for(ind.name in ind.names){
  res.mean <- paste0(pth2valid, ind.name,'.ReSu.only.rep_test_validation.csv') %>%
    read_csv(show_col_types = F) %>% select(-rep) %>%
    group_by(.validation) %>%
    summarise_all(mean) %>%
    mutate_if(is.numeric, round, 2)
  names(res.mean)[2:length(res.mean)] <- paste0('mean.', names(res.mean)[2:length(res.mean)])
  res.sd <- paste0(pth2valid, ind.name,'.ReSu.only.rep_test_validation.csv') %>%
    read_csv(show_col_types = F) %>% select(-rep) %>%
    group_by(.validation) %>%
    summarise_all(sd) %>%
    mutate_if(is.numeric, round, 4)
  names(res.sd)[2:length(res.sd)] <- paste0('sd.', names(res.sd)[2:length(res.sd)])
  res[[ind.name]] <- left_join(res.mean, res.sd, '.validation')
}
res <- res %>%
  bind_rows(.id='EIV') %>%
  arrange(.validation) %>%
  filter(.validation!='external_static') %>%
  mutate(RMSE = paste0(mean.rmse, ' (SD: ',sd.rmse,')')) %>%
  mutate(rsq = paste0(mean.rsq, ' (SD: ',sd.rsq,')')) %>%
  mutate(cor = paste0(mean.cor, ' (SD: ',sd.cor,')')) %>%
  select(-contains('mean.'), -contains('sd.'))

# export results in a table
res %>%
  flextable::flextable() %>% 
  flextable::save_as_docx( path = "./code/data/validation/ReSu.only.rep_test_validation.docx")

# quit
quit(save = 'no')