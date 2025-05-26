################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Tune and train Random Forests using EUNIS-habitat level 2

################################################################################

#### 1. Data preparation and model set up ####

# select response variable name ('Moisture by default')
ind.name <- 'EIV_M' # use `EIV_N` for nutrients, `EIV_T` for temperature, etc.

# folder to store model results
pth2export <- './models/ESy2models/'

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


# extract ESy2 and discard unwanted plots based on EUNIS-ESy level-2 habitat classification
dat_split$ESy2 <- # sapply extract_habitat_lev2() function to extract EUNIS-ESy level-2 habitat code
  unname(sapply(dat_split$ESy, extract_habitat_lev2)) 
dat_split <- dat_split %>%
  filter(!is.na(ESy2)) %>% # remove all habitats without EUNIS-ESY level 2
  filter(ESy2 != 'N1') %>% # remove 'Coastal dunes and sandy shores'
  filter(ESy2 != 'Q3')     # remove 'Palsa mires' - they are nearly all out of study region

table(dat_split$ESy2) # list of habitats we will focus on

# rename target variables (treshold of species for plot inclusion and taret CMeiv variable)
names(dat_split)[which(names(dat_split) == paste0('n.', ind.name))]  <- 'treshold'
names(dat_split)[which(names(dat_split) == paste0('cm.', ind.name))] <- 'eiv'

# apply filters accordingly to the raw data
treshold.of.EIVE.species = 0.8 # set minimum proportion of species with available EIV to include in the analyses
dat_split <- dat_split %>%
  filter(treshold >= treshold.of.EIVE.species)

# final selection of variables needed for modeling
dat_split <- dat_split %>%
  select(plot_id, eiv, x, y, elev, year, ESy2, n, plot_size)

# #work on a random subset?
# dat_split <- dat_split %>% sample_n(5000)

# split the data
dat_split <- dat_split %>%
  initial_split(prop = 4 / 5,  # We will use 80% of the data for training, 20% for testing
                strata = ESy2  # Ensure ESy2 habitats will be sampled in both training and testing set
                )

# subset training set
dat_train <- training(dat_split)

# get CV folds
set.seed(234)
cv_random_folds <-
  vfold_cv(dat_train,
           v = 10,
           repeats = 1,
           strata = ESy2 # Ensure ESy2 habitats will be sampled in both training and testing set
           )

# define recipe
rec <- recipe(eiv ~ x + y + elev + year + ESy2 + n + plot_size, data = dat_train) 

# define model
spec <- rand_forest(trees = 500,
                    min_n = tune(),
                    mtry = tune()) %>%
  set_mode('regression') %>%
  set_engine('ranger', importance = 'impurity', seed = 1975)

# define workflow
wflow <- workflow() %>%
  add_model(spec) %>%
  add_recipe(rec)

# define tune grid
grd_tune <- expand.grid(
  min_n = c(2, 5, 10, 20),
  mtry = c(2:7)) %>%
  as_tibble()

#### 2. Perform tuning ####
# set up paralleling process
plan(multisession, workers = 10)

set.seed(345)
st = Sys.time()
tune_res <- tune_grid(
  object = wflow,
  resamples = cv_random_folds,
  grid = grd_tune,
  # control = control_grid(verbose = T),
  metrics = metric_set(rmse, rsq)
)
print(Sys.time() - st)

# export tuning results
tune_res %>%
  write_rds(paste0(pth2export, 'RF.tune_res_', ind.name, '.rds'),
            compress = 'gz')

# display the best performing sets of parameters
show_best(tune_res, metric = 'rmse')

#### 3. Fit models with best hyperparameters ####

# Select the best performing sets of parameters
tune_best <- select_best(tune_res, metric = 'rmse')
tune_best

# Finalize workflow, fit on training data and test on testing data
lfit <- wflow %>%
  finalize_workflow(tune_best) %>%
  last_fit(dat_split)

# Export last fit
lfit %>%
  write_rds(paste0(pth2export, 'RF.last_fit_', ind.name, '.rds'),
            compress = 'gz')

# Collect eval. metrics last fit
lfit %>%
  collect_metrics()

# stop multisession
plan(sequential)

# quit
quit(save = 'no')