################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Anonymize, clean, and zip EVA and ReSurveyEU plot records 
# This is done to share data in the repository
# N.B. The data available in the repository already went through this process!
# DO NOT RUN

################################################################################

library(tidyverse)

# Prepare EVA data ####

EVA <- './data/EVA.csv' %>%
  read_csv(show_col_types = F) %>%
  select(-contains('sd.')) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(elev = round(elev, 0)) %>%
  arrange(desc(x), desc(y)) %>%
  mutate(plot_id = 1:nrow(.)) # anonymize plot id

# anonymize datasets
EVA_datasets <- EVA %>% 
  select(dataset) %>%
  unique() %>%
  mutate(dataset_abbr =
           abbreviate(stringi::stri_trans_general(dataset, 'Latin-ASCII'), 4, named = F) %>%
            str_replace_all('_','') %>% str_replace_all('-','') %>% str_replace_all('\\.','') %>%
            tolower()
         ) 
EVA_datasets$dataset %>% duplicated() %>% table
EVA_datasets$dataset_abbr %>% duplicated() %>% table

EVA <- EVA %>% 
  left_join(EVA_datasets, 'dataset') %>%
  select(1,2, dataset_abbr, everything()) %>%
  select(-dataset) %>%
  rename(dataset = dataset_abbr)

# Prepare ReSurvey data ####

ReSuEU <- './data/ReSurveyEU_clean.csv' %>%
  read_csv(show_col_types = F) %>%
  select(-contains('sd.'), -ReSur_time, -ReSur_obs, -ReSur_plot, -ReSur_site) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(elev = round(elev, 0)) %>%
  arrange(desc(x), desc(y)) %>%
  mutate(plot_id = seq((max(EVA$plot_id)+1), (max(EVA$plot_id))+nrow(.), 1)) # anonymize plot id, continuing EVA series

# anonymize datasets
ReSuEU_datasets <- ReSuEU %>% 
  select(dataset) %>%
  unique() %>%
  mutate(dataset_abbr =
           abbreviate(stringi::stri_trans_general(dataset, 'Latin-ASCII'), 6, named = F) %>%
           str_replace_all('_','') %>% str_replace_all('-','') %>% str_replace_all('\\.','') %>%
           tolower()
  ) 
ReSuEU_datasets$dataset %>% duplicated() %>% table
ReSuEU_datasets$dataset_abbr %>% duplicated() %>% table
table(ReSuEU_datasets$dataset_abbr %in% EVA$dataset) # check no overlap in EVA

ReSuEU <- ReSuEU %>% 
  left_join(ReSuEU_datasets, 'dataset') %>%
  select(1:4, dataset_abbr, everything()) %>%
  select(-dataset) %>%
  rename(dataset = dataset_abbr)

# anonymize Resurvey id
ReSuEU_resurv_id <- ReSuEU %>% 
  select(resurv_id) %>%
  unique() %>%
  mutate(resurv_id2=1:nrow(.))

ReSuEU <- ReSuEU %>% 
  left_join(ReSuEU_resurv_id, 'resurv_id') %>%
  select(1, resurv_id2, everything()) %>%
  select(-resurv_id) %>%
  rename(resurv_id = resurv_id2) %>%
  # arrange everything by id and year
  arrange(resurv_id, year)

# reorder a bit
ReSuEU <- select(ReSuEU, database:n, x_mean, y_mean, max_dist_m, everything())

#Export `EVA` and `ReSurveyEU_clean` results####
write_csv(EVA, './data/EVA.csv.xz')
write_csv(ReSuEU, './data/ReSurveyEU_clean.csv.xz')


# Anonymize, clean, and zip EVA and ReSurveyEU subset for CWM ####
EVA_ReSu_CWM <- './data/EVA_ReSu_CWM.csv' %>%
  read_csv(show_col_types = F) %>%
  select(-contains('sd.')) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(elev = round(elev, 0)) %>%
  arrange(desc(x), desc(y)) %>%
  mutate(plot_id = 1:nrow(.)) # anonymize plot id

# anonymize datasets
EVA_ReSu_CWM_datasets <- EVA_ReSu_CWM %>% 
  select(dataset) %>%
  unique() %>%
  mutate(dataset_abbr =
           abbreviate(stringi::stri_trans_general(dataset, 'Latin-ASCII'), 7, named = F) %>%
           str_replace_all('_','') %>% str_replace_all('-','') %>% str_replace_all('\\.','') %>%
           tolower()
  ) 
EVA_ReSu_CWM_datasets$dataset %>% duplicated() %>% table
EVA_ReSu_CWM_datasets$dataset_abbr %>% duplicated() %>% table

EVA_ReSu_CWM <- EVA_ReSu_CWM %>% 
  left_join(EVA_ReSu_CWM_datasets, 'dataset') %>%
  select(1,2, dataset_abbr, everything()) %>%
  select(-dataset) %>%
  rename(dataset = dataset_abbr)

# Export
write_csv(EVA_ReSu_CWM, './data/EVA_ReSu_CWM.csv.xz')

# Anonymize, clean, and zip EVA and ReSurveyEU subset for CM calculated without trees ####
EVA_ReSu_NOTREES <- './data/EVA_ReSu_NOTREES.csv' %>%
  read_csv(show_col_types = F) %>%
  select(-contains('sd.')) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(elev = round(elev, 0)) %>%
  arrange(desc(x), desc(y)) %>%
  mutate(plot_id = 1:nrow(.)) # anonymize plot id

# anonymize datasets
EVA_ReSu_NOTREES_datasets <- EVA_ReSu_NOTREES %>% 
  select(dataset) %>%
  unique() %>%
  mutate(dataset_abbr =
           abbreviate(stringi::stri_trans_general(dataset, 'Latin-ASCII'), 7, named = F) %>%
           str_replace_all('_','') %>% str_replace_all('-','') %>% str_replace_all('\\.','') %>%
           tolower()
  ) 
EVA_ReSu_NOTREES_datasets$dataset %>% duplicated() %>% table
EVA_ReSu_NOTREES_datasets$dataset_abbr %>% duplicated() %>% table

EVA_ReSu_NOTREES <- EVA_ReSu_NOTREES %>% 
  left_join(EVA_ReSu_NOTREES_datasets, 'dataset') %>%
  select(1,2, dataset_abbr, everything()) %>%
  select(-dataset) %>%
  rename(dataset = dataset_abbr)

# Export
write_csv(EVA_ReSu_NOTREES, './data/EVA_ReSu_NOTREES.csv.xz')