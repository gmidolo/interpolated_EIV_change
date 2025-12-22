################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 22.12.2025
################################################################################

# Description: Interpolate predictions and compare CWM vs. CM trends

################################################################################

#1. Interpolation ####

# load packages
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(tidymodels)
    library(matrixStats)
  }
)

# dir where models are stored
pth2models <- './models/'

# dir where pdp data should be stored
pth2pdp <- './preds/pdp/CWM.vs.CM/'

# dir where figures are stored
pth2fig <- './fig/'

# standardize plot size to median by habitat?
standardize_plot.size = TRUE

# define range of predictions
pred_years <- seq(1960, 2020, 1) 

# load plot data for interpolation
set.seed(123)
dat.i <- './data/EVA_ReSu_CWM.csv.xz' %>% # load subset data containing CWM
  read_csv(show_col_types = F) %>%
  split(., .$database)
dat.i$ReSurveyEU <- dat.i$ReSurveyEU %>% # format ReSurveyEurope data
  arrange(resurv_id, year) %>%
  group_by(resurv_id) %>%
  slice_sample(n=1) %>%
  ungroup() 
dat.i <- dat.i %>%
  bind_rows()

# calculate habitat plot size medians
if (standardize_plot.size) {
  hab_plot.size <- dat.i %>%
    group_by(dataset, habitat) %>%
    summarise(plot_size_median_in_db = median(plot_size),n_in_db=n()) %>%
    ungroup() %>%
    group_by(habitat) %>%
    summarise(
      plot_size_median = weightedMedian(x=plot_size_median_in_db, w=n_in_db)
    )
}

# interpolate over study period
dat.i <- dat.i %>%
  filter(year >= min(pred_years) & year <= max(pred_years))

# set up variable names
ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))
wei.names = c('CM','CWM')

# interpolation loop
for (wei in wei.names) {
  for(ind.name in ind.names$eiv_name_raw){
    
    pred_res <- list()
    
    st <- Sys.time()
    
    print(paste0('Start interpolation for ', ind.name, ' - ', wei))
    
    #prepare data 
    dat <- dat.i

    # rename
    names(dat)[which(names(dat) == paste0('n.', ind.name))] <- 'threshold'
    names(dat)[which(names(dat) == paste0(tolower(wei), '.', ind.name))] <- 'eiv'
    
    # apply filters accordingly to the raw data
    threshold.of.EIVE.species = 0.8 # set minimum proportion of species with available EIV to include in the analyses
    dat <- dat %>%
      filter(threshold >= threshold.of.EIVE.species)
    
    # final selection of variables needed for modeling
    dat <- dat %>%
      select(plot_id, eiv, x, y, elev, year, habitat, n, plot_size)
    
    # Set output file name to export
    output_file <- paste0(pth2pdp, wei, '_', ind.name, '_habitat_pdp_trends_alltrees.csv')
    
    # Load model and extract workflow
    model_path <- paste0(pth2models, wei,'_RF.last_fit_', ind.name, '.rds')
    m <- read_rds(model_path)
    m_wf <- extract_workflow(m)
    
    for (i in pred_years) {
      
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
    pred_res_df <- pred_res_df %>% 
      mutate(weighting = wei, .before = 1)
    
    # Export results
    write_csv(pred_res_df, output_file)
    
    print(Sys.time() - st)
  }
}


#2. Plot ####

# load pdp data
dat.i <- list()
for (ind.name in ind.names$eiv_name_raw) {
  cm <- paste0(pth2pdp, 'CM', '_', ind.name, '_habitat_pdp_trends_alltrees.csv') %>%
    read_csv(show_col_types = F)
  cwm <- paste0(pth2pdp, 'CWM', '_', ind.name, '_habitat_pdp_trends_alltrees.csv') %>%
    read_csv(show_col_types = F)
  dat.i[[ind.name]] <- bind_rows(cm, cwm)
}
dat.i <- dat.i %>% 
  bind_rows(.id = 'eiv_name_raw') %>% 
  left_join(ind.names, by = join_by(eiv_name_raw))

# Isolate means from the 60s
dat.i.60 <- dat.i %>% 
  filter(year == 1960) %>% 
  rename(mean_60 = mean) %>% 
  select(eiv_name, weighting, habitat, mean_60) %>% 
  unique()
dat.i <- dat.i %>%
  left_join(dat.i.60, by = join_by(weighting, habitat, eiv_name)) %>%
  mutate(stndrd_mean = mean - mean_60)

dat.i$eiv_name <- factor(dat.i$eiv_name, ind.names$eiv_name)

# define limits on y axis
y_breaks <- seq(-0.4, 0.4, 0.2)
y_limits_centrd <- range(y_breaks)
# define y axis lables 
y_left <- y_breaks[1:floor(length(y_breaks)/2)]
y_right <- y_breaks[(ceiling(length(y_breaks)/2)+1):length(y_breaks)]
y_right <- paste0('+', y_right)
y_labels <- c(y_left, '0', y_right)

pi <- ggplot(dat.i, aes(year, stndrd_mean, col = weighting)) +
  # plot zero change line
  geom_hline(yintercept = 0, lty=3, color='grey30', linewidth=.3) +
  # plot pdp curves
  geom_line(alpha = .8) +
  # details
  labs(y=expression(paste('CM/CWM', ''[EIV], ' change since 1960')), x='Year', color='Metric') +
  scale_colour_manual(values = c('#1038af','#ff4949')) +
  theme_bw() +
  facet_grid(habitat~eiv_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # set appropriate scale 
  scale_y_continuous(
    limits = c(-.42, .42),
    labels = y_labels,
    breaks = y_breaks
  )
pi

# export figure
# ggsave(paste0(pth2fig, 'CMvsCWM.trends.pdf'), pi, width = 6.8, height = 5)
ggsave(paste0(pth2fig, 'CMvsCWM.trends.png'), pi, width = 6.8, height = 5, dpi = 600)

# quit
quit(save = 'no')