################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 21.12.2025
################################################################################

# Description: Tune and train Random Forests for both CWM and CM on plot subset

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

# CWM (= TRUE) vs. CM (= FALSE)
weighted.values <- c(TRUE, FALSE)

# set number of workers (one for each EIV and weighting combination)
no.total.workers <- 10

# plan workers (outer parallelization only)
plan(multisession, workers = no.total.workers)

# get model grid
model_grid <- expand_grid(
  ind.name = ind.name.values,
  weighted = weighted.values
)

# plan workers
plan(multisession, workers = no.total.workers)

# assign each combination to workers
model_grid_rows <- split(model_grid, seq_len(nrow(model_grid)))

# map in parallel
future_map(model_grid_rows, function(params) {
  
  # 1. Data preparation and model set up 
  
  # select response variable name
  ind.name <- params$ind.name
  
  # model either weighted CWM or unweighted CM
  weighted <- params$weighted
  test.name <- ifelse(weighted, 'CWM', 'CM')
  
  # export folder
  pth2export <- './models/'
  
  # load only packages inside the worker
  suppressPackageStartupMessages({
    library(tidyverse)
    library(tidymodels)
    library(future)
  })
  
  # prepare data
  set.seed(123)
  
  dat_split <- './data/EVA_ReSu_CWM.csv.xz' %>%
    read_csv(show_col_types = FALSE) %>%
    split(., .$database)
  
  dat_split$ReSurveyEU <- dat_split$ReSurveyEU %>%
    arrange(resurv_id, year) %>%
    group_by(resurv_id) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  dat_split <- bind_rows(dat_split)
  
  # rename vars
  names(dat_split)[which(names(dat_split) == paste0('n.', ind.name))] <- 'threshold'
  
  if (weighted) {
    names(dat_split)[which(names(dat_split) == paste0('cwm.', ind.name))] <- 'eiv'
  } else {
    names(dat_split)[which(names(dat_split) == paste0('cm.', ind.name))] <- 'eiv'
  }
  
  # filter
  threshold.of.EIVE.species = 0.8
  dat_split <- dat_split %>% filter(threshold >= threshold.of.EIVE.species)
  
  # final vars
  dat_split <- dat_split %>%
    select(plot_id, eiv, x, y, elev, year, habitat, n, plot_size)
  
  # split
  dat_split <- dat_split %>%
    initial_split(prop = 4/5, strata = eiv)
  dat_train <- training(dat_split)
  
  # CV folds
  set.seed(234)
  cv_random_folds <- vfold_cv(dat_train, v = 10, repeats = 1, strata = eiv)
  
  # recipe
  rec <- recipe(eiv ~ x + y + elev + year + habitat + n + plot_size, data = dat_train) %>%
    step_string2factor(habitat, levels = c('Forest', 'Grassland', 'Scrub', 'Wetland'))
  
  # model
  spec <- rand_forest(
    trees = 500,
    min_n = tune(),
    mtry = tune()
  ) %>%
    set_mode('regression') %>%
    set_engine('ranger', importance = 'impurity', seed = 1975)
  
  # workflow
  wflow <- workflow() %>%
    add_model(spec) %>%
    add_recipe(rec)
  
  # tuning grid
  grd_tune <- expand.grid(
    # use a smaller subset of parameter grids
    min_n = c(2,5), 
    mtry = c(3:7)
  ) %>% as_tibble()
  
  # 2. Tuning
  set.seed(345)
  st <- Sys.time()
  tune_res <- tune_grid(
    object = wflow,
    resamples = cv_random_folds,
    grid = grd_tune,
    metrics = metric_set(rmse, rsq)
  )
  print(Sys.time() - st)
  
  # # export tuning
  # tune_res %>%
  #   write_rds(paste0(pth2export, test.name, '_RF.tune_res_', ind.name, '.rds'), compress = 'gz')
  
  show_best(tune_res, metric = 'rmse')
  
  # 3. Fit best model
  tune_best <- select_best(tune_res, metric = 'rmse')
  
  lfit <- wflow %>%
    finalize_workflow(tune_best) %>%
    last_fit(dat_split)
  
  lfit %>%
    write_rds(paste0(pth2export, test.name, '_RF.last_fit_', ind.name, '.rds'), compress = 'gz')
  
  print(paste('Fit completed for:', ind.name, '- Test:', test.name))
  
})

# stop multisession
plan(sequential)

# quit
quit(save='no')