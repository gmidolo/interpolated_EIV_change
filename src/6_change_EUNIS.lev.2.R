################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Interpolate EIV changes across EUNIS-ESy level 2 habitats

################################################################################

#### 1. Set up data ####

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(flextable)
  library(scales)
  library(matrixStats)
})

# standardize plot size to median by habitat?
standardize_plot.size = TRUE

# set indicator names
ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

# source help functions 
source('./src/0_helpfunctions.R')

# folder where to store figures and diagnostic output
pth2fig <- './fig/'

# prepare data for modeling and split train and test dataset
set.seed(123)
dat.initial <-
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope( # Load ReSurveyEU (static, using one random point in the survey)
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
    )[['traintest_data']]  
  )

# standardize plot size to median by habitat while predicting EIV values?
standardize_plot.size = TRUE

# n of plots with complete information
dat.initial %>%
 filter(n.EIV_M >= .8 | n.EIV_N >= .8 | n.EIV_T >= .8 | n.EIV_L >= .8 | n.EIV_R >= .8 ) %>%
 nrow()

# calc year range for interpolation for each habitat
dat.year.ESy2 <- dat.initial
treshold.of.EIVE.species = 0.8
dat.year.ESy2 <- dat.initial %>%
  filter(n.EIV_M >= treshold.of.EIVE.species) %>%
  filter(n.EIV_N >= treshold.of.EIVE.species) %>%
  filter(n.EIV_T >= treshold.of.EIVE.species) %>%
  filter(n.EIV_L >= treshold.of.EIVE.species) %>%
  filter(n.EIV_R >= treshold.of.EIVE.species)
dat.year.ESy2$ESy2 <-
  unname(sapply(dat.year.ESy2$ESy, extract_habitat_lev2))
dat.year.ESy2 <- dat.year.ESy2
dat.year.ESy2 <- dat.year.ESy2 %>%
  group_by(ESy2) %>%
  summarise(max_yr = round(quantile(year, .95)),
            min_yr = round(quantile(year, .05)))
head(dat.year.ESy2)

# Interpolate changes using 2020 vs 1960? If FALSE get data back to 2020 - 1960 period
interpolate_to_2020.1960 = FALSE
if(interpolate_to_2020.1960) {
  dat.year.ESy2$max_yr = 2020
  dat.year.ESy2$min_yr = 1960
}

# function for 0-overlap test (calc. significant vs. non-significant CI)
CI_overlap_test <- function(x, y) {
  zero_check <- (x == 0 & y == 0)
  sign_check <- (x > 0 & y > 0) | (x < 0 & y < 0)
  res <- zero_check | sign_check
  return(res)
}

#### 2. Predict ####
# obtain average percentage changes per EUNIS-ESy habitat lev.2 and 95%CI
for(ind.name in ind.names$eiv_name_raw) {
  print(paste0('Start interpolation for ', ind.name))
  
  st = Sys.time()
  dat = dat.initial
  
  # assign names to focal EIV variables
  names(dat)[which(names(dat) == paste0('n.', ind.name))] <- 'treshold'
  names(dat)[which(names(dat) == paste0('cm.', ind.name))] <- 'eiv'
  
  # apply filters accordingly to the raw data
  treshold.of.EIVE.species = 0.8
  dat <- dat %>%
    filter(treshold >= treshold.of.EIVE.species)
  
  # extract ESy2 and discard unwanted plots based on EUNIS-ESy level-2 habitat classification
  dat$ESy2 <- unname(sapply(dat$ESy, extract_habitat_lev2))
  dat <- dat %>%
    filter(!is.na(ESy2)) %>%
    filter(ESy2 != 'N1') %>%
    filter(ESy2 != 'Q3')
  
  # final variable selection
  dat <- dat %>%
    select(plot_id, dataset, eiv, x, y, elev, year, ESy2, n, plot_size)
  
  # Calculate habitat plot size medians
  if (standardize_plot.size) {
    hab_plot.size <- dat %>%
      group_by(dataset, ESy2) %>%
      summarise(plot_size_median_in_db = median(plot_size),
                n_in_db = n()) %>%
      ungroup() %>%
      group_by(ESy2) %>%
      summarise(plot_size_median = weightedMedian(x = plot_size_median_in_db, w =
                                                    n_in_db))
  }
  
  # load model fit (ESy2-only models)
  pth_model <-
    list.files(
      './models/ESy2models/',
      pattern = paste0('RF.last_fit_', ind.name),
      full.names = T
    )
  m <- pth_model %>% read_rds()
  
  # extract workflow for predictions
  m_wf <- extract_workflow(m)
  
  # predict EIVs
  pred_dat <- dat
  for (i in c('min_yr', 'max_yr')) {
    cat('predicting over ', i, ' - ')
    pdi <- pred_dat
    pdi <- pdi %>%
      select(-year) %>%
      left_join(dat.year.ESy2[, names(dat.year.ESy2) %in% c('ESy2', i)], by =
                  'ESy2')
    names(pdi)[names(pdi) == i] <- 'year'
    if (standardize_plot.size) {
      # standardize plot size to median by habitat?
      pdi <- pdi %>%
        select(-plot_size) %>%
        left_join(hab_plot.size, 'ESy2') %>%
        rename(plot_size = plot_size_median)
    }
    pred_dat[paste0('eiv_pred_', i)] <- predict(object = m_wf, new_data = pdi) %>% pull(.pred)
  }
  
  # predict changes
  pred_focalchange <- pred_dat %>%
    mutate(eiv_abs.change = eiv_pred_max_yr - eiv_pred_min_yr) %>%
    select(-contains('_pred_'))
  
  # summarize stats
  res.i <- pred_focalchange %>%
    group_by(ESy2) %>%
    summarise(
      mean = mean(eiv_abs.change),
      sd = sd(eiv_abs.change),
      n = n()
    ) %>%
    ungroup() %>%
    mutate(se = sd / sqrt(n)) %>%
    mutate(lci = mean - (1.96 * se),
           uci = mean + (1.96 * se))
  res.i$significance <- CI_overlap_test(res.i$lci, res.i$uci)
  
  # finalize table
  res.i <- res.i %>%
    mutate_if(is.numeric, round, 2) %>%
    mutate(EIV = ind.names[ind.names$eiv_name_raw == ind.name, ]$eiv_name, .before =
             ESy2)
  
  # store results in the list
  res[[ind.name]] <- res.i
  
}

#### 3. Format final table and export ####
getESy1 <- \(ESy2.code){
  str=substr(ESy2.code,1,1)
  if(str == 'T'){return('Forest')}
  if(str == 'R'){return('Grassland')}
  if(str == 'S'){return('Scrub')}
  if(str == 'Q'){return('Wetland')}
}
res1 <- res %>%
 map(\(x){
  x %>%
  mutate(ci = uci - mean) %>%
  mutate(sign = ifelse(.$significance, '', 'n.s.')) %>%
  mutate(res.sum = paste0(mean, '±(', round(ci, 2),')', sign)) %>%
  select(ESy2, EIV, res.sum) %>%
  spread(EIV, res.sum) %>%
  mutate(`Habitat type` = unname(sapply(.$ESy2, getESy1)), .before=ESy2) %>%
  arrange(`Habitat type`)
 })

res_n <- bind_rows(res) %>% 
  group_by(ESy2) %>% summarise(n_median = median(n)) %>% ungroup() %>%
  mutate(`Habitat type` = unname(sapply(.$ESy2, getESy1)), .before=ESy2) %>%
  arrange(`Habitat type`) 

res_fin <- res_n %>%
           left_join(res1$EIV_M) %>%
           left_join(res1$EIV_N) %>%
           left_join(res1$EIV_T) %>%
           left_join(res1$EIV_L) %>%
           left_join(res1$EIV_R)

pretty_years <- dat.year.ESy2 %>%
 mutate(Years = paste0(min_yr,'-',max_yr)) %>%
 select(ESy2, Years)

res_fin <- res_fin %>%
  left_join(pretty_years, 'ESy2') %>%
  select(1:2, Years, everything())

# load pretty ESy2-names
hab.names <- readxl::read_excel('EUNIS_ESy2_habitat.names.xlsx', col_names = F) %>%
  setNames('Habitat name') %>%
  mutate(ESy2 = substr(`Habitat name`, 1, 2)) %>%
  mutate(`Habitat name`=substr(`Habitat name`, 4, nchar(`Habitat name`)))

res_fin <- res_fin %>%
  left_join(hab.names) %>%
  select(1:2, `Habitat name`, everything())

# export
res_fin %>%
 select(-ESy2) %>%
 rename(`no. plots`=n_median) %>%
 flextable::flextable() %>% 
 flextable::save_as_docx(path = paste0(pth2fig, 'habitat_ESy2_means.docx'))

# quit
quit(save = 'no')