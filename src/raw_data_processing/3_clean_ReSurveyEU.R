#--- Clean ReSurveyEU ---#

# Here, we 1) filter surveys with at least two plots and 2) filter resurvey with at least two plots sampled in different years

# Load packages
suppressPackageStartupMessages(
  {library(tidyverse); library(sf)}
)

username <- str_split(getwd(), '/')[[1]][3]
setwd(paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU_EIV/'))

date.update = '01.11.2025'

resurvey_clean <- paste0('./code/data/input/ReSurveyEU','_',date.update,'.csv') %>%
                   read_csv(show_col_types = F)   %>%
                  # define putative re-survey id (group plots within same survey):
                  group_by(ReSur_type, ReSur_site, ReSur_plot) %>%
                  mutate(resurv_id = cur_group_id(), .before = plot_id) %>% 
                  ungroup() %>%
                  # filter surveys with at least two plots:
                  group_by(resurv_id) %>%
                  filter(n()>=2) %>%
                  # filter resurvey with at least two samplings in different years
                  filter(length(unique(year)) >= 2) %>%
                  arrange(resurv_id, year) %>%
                  ungroup()

# Synthetize coordinates
resurv_id_list <- resurvey_clean %>% split(., .$resurv_id)
coords_synthesis <- list()
for(i in names(resurv_id_list)){
 crdi <- resurv_id_list[[i]][, c('x','y')]
 summ <- summarise_all(crdi, mean) %>%
         setNames(c('x_mean','y_mean'))
coords_synthesis[[i]] <- cbind(summ, data.frame(max_dist_m = max(dist(crdi, diag=T))))
}
coords_synthesis <- bind_rows(coords_synthesis, .id='resurv_id')
plot(coords_synthesis$x_mean, coords_synthesis$y_mean)

resurvey_clean <- resurvey_clean %>%
                  left_join(coords_synthesis %>% mutate(resurv_id=as.numeric(resurv_id)), 'resurv_id')

# Remove plots with uncertainty in plots relocation higher than 100 m
prop.table(table(coords_synthesis$max_dist_m > 100))
nrow(resurvey_clean)
resurvey_clean <- resurvey_clean %>%
 filter(max_dist_m <= 100)
nrow(resurvey_clean) # total number of plots in the ReSurveyEU database
table(resurvey_clean$ReSur_type) # Resurvey type

resurvey_clean$resurv_id %>% unique() %>% length() # number of resurveys (time series)
resurvey_clean$plot_id %>% unique() %>% length()   # number of plots

# Export
resurvey_clean %>%
 write_csv(
   paste0('./code/data/input/ReSurveyEU_clean','_',date.update,'.csv')
 )