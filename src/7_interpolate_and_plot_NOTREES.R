################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 22.12.2025
################################################################################

# Description: Interpolate predictions and compare analysis without trees/shrubs

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
pth2pdp <- './preds/pdp/no.trees/'

# dir where figures are stored
pth2fig <- './fig/'

# standardize plot size to median by habitat?
standardize_plot.size = TRUE

# define range of predictions
pred_years <- seq(1960, 2020, 1) 

# load plot data for interpolation
set.seed(123)
dat.i <- './data/EVA_ReSu_NOTREES.csv.xz' %>%
  read_csv(show_col_types = F) %>%
  split(., .$database)
dat.i$ReSurveyEU <- dat.i$ReSurveyEU %>%
  arrange(resurv_id, year) %>%
  group_by(resurv_id) %>%
  slice_sample(n=1) %>%
  ungroup() 
dat.i <- dat.i %>% bind_rows()

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

# set up response variable names
ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))
comparison.names <- c('noTREES','all')

# focus on forests and scrub vegetation only
dat.i <- dat.i %>% 
  semi_join(data.frame(habitat = c('Forest', 'Scrub')), by = 'habitat')

# prediction loop
res <- list()
i <- 1
for (comparison.name in comparison.names) {
  for (ind.name in ind.names$eiv_name_raw) {
    
    pred_res <- list()
    st <- Sys.time()
    
    print(paste0('Start interpolation for ', ind.name, ' - ', comparison.name))
    
    # prepare data 
    dat <- dat.i
    
    # rename cols
    names(dat)[which(names(dat) == paste0('n.', ind.name))]  <- 'threshold'
    names(dat)[which(names(dat) == paste0('cm.', ind.name))] <- 'eiv'
    
    # filter
    threshold.of.EIVE.species = 0.8
    dat <- dat %>% 
      filter(threshold >= threshold.of.EIVE.species) %>%
      select(plot_id, eiv, x, y, elev, year, habitat, n, plot_size)
    
    # load appropriate model
    model_path <- if (comparison.name == 'noTREES') {
      paste0(pth2models, 'NOTREES_RF.last_fit_', ind.name, '.rds')
    } else {
      paste0(pth2models, 'RF.last_fit_', ind.name, '.rds')
    }
    m_wf <- read_rds(model_path) %>% 
      extract_workflow()
    
    # predictions over years 
    for (i in pred_years) {
      pdi <- dat %>% mutate(year = i)
      
      if (standardize_plot.size) {
        pdi <- pdi %>% 
          select(-plot_size) %>% 
          left_join(hab_plot.size, 'habitat') %>%
          rename(plot_size = plot_size_median)
      }
      
      praw <- predict(m_wf, new_data = pdi, type = 'raw', opts = list(predict.all = TRUE))
      praw <- as_tibble(praw)
      
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
          n = n(),
          sd = sd(.pred),
          mean = mean(.pred),
          lower_ci_95 = mean - (1.96 * (sd/sqrt(n))),
          upper_ci_95 = mean + (1.96 * (sd/sqrt(n))),
          .groups = 'drop'
        )
    }
    
    # bind yearly predictions
    pred_res_df <- bind_rows(pred_res, .id = 'year') %>%
      mutate(
        comparison = comparison.name,
        eiv_name = ind.names %>% filter(eiv_name_raw == ind.name) %>% pull(eiv_name),
        .before = 1
      )
    
    # store in list
    res[[i]] <- pred_res_df
    i <- i + 1
    
    print(Sys.time() - st)
  }
}

# combine all results into one dataframe
final_results <- bind_rows(res)

# export combined output
output_file <- paste0(pth2pdp, 'NOTREES_pdp_trends_combined.csv')
write_csv(final_results, output_file)


#2. Plot ####

# load pdp data 
d <- paste0(pth2pdp, 'NOTREES_pdp_trends_combined.csv') %>%
  read_csv(show_col_types = F)

# Isolate means from the 60s
d.60 <- d %>% 
  filter(year == 1960) %>%
  rename(mean_60 = mean) %>% 
  select(eiv_name, comparison, habitat, mean_60) %>% 
  unique()
d <- d %>%
  left_join(d.60, by = join_by(comparison, eiv_name, habitat)) %>%
  mutate(stndrd_mean = mean - mean_60)

d$eiv_name <- factor(d$eiv_name, ind.names$eiv_name)

# define limits on y axis
y_breaks <- seq(-0.2, 0.4, 0.2)
y_limits_centrd <- range(y_breaks)
# define y axis lables 
y_left <- y_breaks[1:floor(length(y_breaks)/2)]
y_right <- y_breaks[(ceiling(length(y_breaks)/2)+1):length(y_breaks)]
y_right <- paste0('+', y_right)
y_labels <- c(y_left, '0', y_right)

pi <- ggplot(d, aes(year, stndrd_mean, col = comparison)) +
  # plot zero change line
  geom_hline(yintercept = 0, lty=3, color='grey30', linewidth=.3) +
  # plot pdp curves
  geom_line(alpha = .8) +
  # details
  labs(y=expression(paste('CM', ''[EIV], ' change since 1960')), x='Year', color='Metric') +
  scale_color_manual(
    labels = c('all' = 'All species', 'noTREES' = 'Tree species excluded'),
    values = c('all' = '#1038af', 'noTREES' = '#ff4949') 
  ) +
  theme_bw() +
  facet_grid(habitat~eiv_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # set appropriate scale 
  scale_y_continuous(
    limits = y_limits_centrd,
    breaks = y_breaks
  )
pi

# export figure
# ggsave(paste0(pathtofig, 'all.vs.NOTREES.trends.pdf'), pi, width = 7.25, height = 3.75)
ggsave(paste0(pth2fig, 'all.vs.NOTREES.trends.png'), pi, width = 7.25, height = 3.75, dpi = 600)

# quit
quit(save = 'no')