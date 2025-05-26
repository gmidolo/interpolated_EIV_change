################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Apply 10-fold CV with 3 repeats with Random Forests

################################################################################

#### 1. Data preparation and model set up ####

# select response variable name ('Moisture by default')
ind.name <- 'EIV_M' # use `EIV_N` for nutrients, `EIV_T` for temperature, etc.

# folder to store model results
pth2export <- './models/'

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(future)
})

# source functions to format Resurvey data
source('./src/0_helpfunctions.R')

# prepare data for modeling and split train and test dataset
set.seed(123)
dat_split <-
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope( # Load ReSurveyEU (static, using one random point in the survey)
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
    )[['traintest_data']]  
  )

# rename target variables (treshold of species for plot inclusion and taret CMeiv variable)
names(dat_split)[which(names(dat_split) == paste0('n.', ind.name))]  <- 'treshold'
names(dat_split)[which(names(dat_split) == paste0('cm.', ind.name))] <- 'eiv'

# apply filters accordingly to the raw data
treshold.of.EIVE.species = 0.8 # set minimum proportion of species with available EIV to include in the analyses
dat_split <- dat_split %>%
  filter(treshold >= treshold.of.EIVE.species)

# final selection of variables needed for modeling
dat_split <- dat_split %>%
  select(plot_id, eiv, x, y, elev, year, habitat, n, plot_size)

# #work on a random subset?
# dat_split <- dat_split %>% sample_n(5000)

# split the data
dat_split <- dat_split %>%
  initial_split(prop = 4 / 5, strata = eiv) # We will use 80% of the data for training, 20% for testing

# subset training set
dat_train <- training(dat_split)


#### 2. Random CV ####

# Load results of the tuning procedure
tune_res <- paste0(pth2export, 'RF.tune_res_', ind.name, '.rds') %>%
  read_rds()

# Select best set of hyperparameters (based on RMSE)
tune_best <- select_best(tune_res, metric = 'rmse')

# Define recipe
rec <-
  recipe(eiv ~ x + y + elev + year + habitat + n + plot_size, data = dat_train) %>%
  step_string2factor(habitat, levels = c('Forest', 'Grassland', 'Scrub', 'Wetland'))

# Define model
spec <- rand_forest(trees = 500,
                    min_n = tune(),
                    mtry = tune()) %>%
  set_mode('regression') %>%
  set_engine('ranger', importance = 'impurity', seed = 1975)

# Define workflow
wflow <- workflow() %>%
  add_model(spec) %>%
  add_recipe(rec) %>%
  finalize_workflow(tune_best)

# Get CV folds
set.seed(124)
cv_random_folds <-
  vfold_cv(dat_train,
           v = 10,
           repeats = 3,
           strata = eiv) # This time we repeat 3 times CV

# Set up paralleling process
plan(multisession, workers = 10)

# Perform CV
set.seed(125)
st = Sys.time()
cv_res <- wflow %>%
  fit_resamples(
    preprocessor = rec,
    resamples = cv_random_folds,
    #control = control_resamples(verbose = T),
    metrics = metric_set(rmse, rsq)
  )
print(Sys.time() - st)

# Export CV raw results
cv_res %>%
  write_rds(paste0(pth2export, 'RF.cv_res_', ind.name, '.rds'), compress = 'gz')

# Export CV metrics
cv_res %>%
  collect_metrics(summarize = F) %>%
  write_csv(paste0(pth2export, 'RF.cv_metrics_', ind.name, '.csv'))

# stop multisession
plan(sequential)

# quit
quit(save = 'no')