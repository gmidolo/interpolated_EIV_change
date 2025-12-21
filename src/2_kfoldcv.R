################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 18.12.2025
################################################################################

# Description: Apply 10-fold CV with 3 repeats with Random Forests

################################################################################

#load packages
suppressPackageStartupMessages({
  library(future)
  library(furrr)
  library(tidyverse)
  library(tidymodels)
})

# indicator names
ind.name.values <- c('EIV_M','EIV_N','EIV_T','EIV_L','EIV_R')

# set number of workers (one for each EIV)
no.total.workers <- length(ind.name.values) 

# plan workers (outer parallelization only)
plan(multisession, workers = no.total.workers)

# map in parallel over indicator names
future_map(ind.name.values, function(ind.name) {
  st = Sys.time()
  
  # 1. data preparation and model setup
  
  # load needed packages inside the worker
  suppressPackageStartupMessages({
    library(tidyverse)
    library(tidymodels)
  })
  
  # export folder
  pth2export <- './models/'
  
  # no.threads in ranger
  no.ranger.workers <- 5
  
  # source function to format Resurvey data
  source('./src/0_helpfunctions.R')
  
  # prepare data for modeling and split train and test dataset
  set.seed(123)
  dat_split <-
    bind_rows(
      # load EVA data
      read_csv('./data/EVA.csv.xz', show_col_types = F),
      # load ReSurveyEU (static, using one random point in the survey)
      format_ReSurveyEurope(
        training_strategy = 'random',
        path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
      )[['traintest_data']]
    )
  
  # rename
  names(dat_split)[which(names(dat_split) == paste0('n.', ind.name))] <- 'treshold'
  names(dat_split)[which(names(dat_split) == paste0('cm.', ind.name))] <- 'eiv'
  
  # apply filters accordingly to the raw data
  treshold.of.EIVE.species = 0.8 # set minimum proportion of species with available EIV to include in the analyses
  dat_split <- dat_split %>%
    filter(treshold >= treshold.of.EIVE.species) 
  
  # final selection of variables needed for modeling
  dat_split <- dat_split %>%
    select(plot_id, eiv, x, y, elev, year, habitat, n, plot_size)
  
  # split the data
  dat_split <- dat_split %>%
    initial_split(prop = 4/5, strata = eiv) # We will use 80% of the data for training, 20% for testing
  
  # subset training set
  dat_train <- training(dat_split) 
  
  # 2. Random CV
  
  # load results of the tuning procedure
  tune_res <- paste0(pth2export,'RF.tune_res_',ind.name,'.rds') %>%
    read_rds()
  
  # select best set of hyperparameters (based on RMSE)
  tune_best <- select_best(tune_res, metric = 'rmse')
  
  # define recipe
  rec <- recipe(eiv ~ x + y + elev + year + habitat + n + plot_size, data = dat_train) %>%
    step_string2factor(habitat, levels = c('Forest', 'Grassland', 'Scrub', 'Wetland'))
  
  # define model 
  spec <- rand_forest(
    trees = 500, 
    min_n = tune(),
    mtry = tune()
  ) %>%
    set_mode('regression') %>%
    set_engine('ranger', importance = 'impurity', seed = 1975, num.threads = no.ranger.workers)
  
  # define workflow
  wflow <- workflow() %>%
    add_model(spec) %>%
    add_recipe(rec) %>%
    finalize_workflow(tune_best)
  
  # get 10 CV folds with 3 repeats 
  set.seed(124)
  cv_random_folds <- vfold_cv(dat_train, v = 10, repeats = 3, strata = eiv) 
  
  # perform CV
  set.seed(125)
  cv_res <- wflow %>%
    fit_resamples(
      preprocessor = rec,
      resamples = cv_random_folds,
      metrics = metric_set(rmse, rsq)
    )
  
  # export CV raw results
  cv_res %>%
    write_rds(paste0(pth2export,'RF.cv_res_', ind.name, '.rds'), compress = 'gz')
  
  # export CV metrics
  cv_res %>%
    collect_metrics(summarize = F) %>%
    write_csv(paste0(pth2export,'RF.cv_metrics_', ind.name, '.csv'))
  
  print(Sys.time()-st)
},
.options = furrr_options(seed = TRUE) # ensure reproducible RNG across futures
)

# stop multisession
plan(sequential)

# quit
quit(save='no')