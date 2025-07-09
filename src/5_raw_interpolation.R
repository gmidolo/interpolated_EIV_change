################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 08.07.2025
################################################################################

# Description: Summarise interpolation stats using raw (all trees) predictions 
#              across main habitat types

################################################################################


# Load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
})

# Set up 
standardize_plot.size <- TRUE
pred_years <- seq(1960, 2020, 1)

# Source function to format Resurvey data
source('./src/0_helpfunctions.R')

# Prepare data for interpolation
set.seed(123)
dat.i <- 
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope(
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz')[['traintest_data']]  # Load ReSurveyEU (static, using one random point in the survey)
  ) 

# Calculate habitat plot size medians
if (standardize_plot.size) {
  hab_plot.size <- dat.i %>%
    group_by(dataset, habitat) %>%
    summarise(plot_size_median_in_db = median(plot_size),n_in_db=n()) %>%
    ungroup() %>%
    group_by(habitat) %>%
    summarise(
      plot_size_median = matrixStats::weightedMedian(x=plot_size_median_in_db, w=n_in_db)
    )
}
hab_plot.size

# Interpolate over study period
dat.i <- dat.i %>%
  filter(year >= min(pred_years) & year <= max(pred_years))

# set up response variable names
ind.names <- c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R')

# Prediction loop
for(ind.name in ind.names){
  
  pred_res <- list()

  st <- Sys.time()
  
  print(paste0('Start interpolation for ', ind.name))

  #### 1. Prepare data ####
  dat = dat.i
  
  # rename
  names(dat)[which(names(dat) == paste0('n.', ind.name))] <- 'treshold'
  names(dat)[which(names(dat) == paste0('cm.', ind.name))] <- 'eiv'

  # apply filters accordingly to the raw data
  treshold.of.EIVE.species = 0.8 # set minimum proportion of species with available EIV to include in the analyses
  dat <- dat %>%
      filter(treshold >= treshold.of.EIVE.species)

  # final selection of variables needed for modeling
  dat <- dat %>%
      select(plot_id, eiv, x, y, elev, year, habitat, n, plot_size)

  # Set output file name to export
  output_file <- paste0('./preds/pdp/', ind.name, '_habitat_pdp_trends_alltrees.csv')

  # Load model and extract workflow
  model_path <- paste0('./models/RF.last_fit_', ind.name, '.rds')
  m <- read_rds(model_path)
  m_wf <- extract_workflow(m)

  for (i in pred_years) {
  
    # message(paste0('Year: ', i))
    
    pdi <- dat %>% mutate(year = i) # change year
    if(standardize_plot.size){ # standardize plot size to median by habitat?
      pdi <- pdi %>% 
        select(-plot_size) %>% 
        left_join(hab_plot.size, 'habitat') %>%
        rename(plot_size = plot_size_median)
      }  
    # get raw ranger predictions
    praw <- predict(m_wf, new_data = pdi, type = 'raw', opts = list(predict.all = TRUE)) 
    praw <- as_tibble(praw) 

    # get stats
    pred_res[[as.character(i)]] <- cbind(pdi, praw) %>%
      group_by(habitat, year) %>%
      summarise(across(starts_with('prediction.'), mean), .groups = 'drop') %>%
      gather('.pred_key', '.pred', starts_with('prediction.')) %>%
      group_by(habitat, year) %>%
      summarise(
        median = median(.pred),
        quantile0.05 = quantile(.pred, .05),
        quantile0.95 = quantile(.pred, .95),
        quantile0.025 = quantile(.pred, .025),
        quantile0.975 = quantile(.pred, .975),
        n=n(),
        sd = sd(.pred),
        mean = mean(.pred),
        lower_ci_95 = mean - (1.96 * (sd/sqrt(n))),
        upper_ci_95 = mean + (1.96 * (sd/sqrt(n))),
      .groups = 'drop'
    )

  }
  
  # Combine predictions across years
  pred_res_df <- bind_rows(pred_res, .id='year')

  # Export results
  write_csv(pred_res_df, output_file)

  print(Sys.time() - st)
}

# Quit
quit(save = 'no')