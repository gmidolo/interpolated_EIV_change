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
  library(cowplot)
  library(ggpubr)
})

# standardize plot size to median by habitat?
standardize_plot.size = TRUE

# set indicator names
ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

# source help functions 
source('./src/0_helpfunctions.R')

# folder where data with predictions (in .csv) were exported
pth2preds <- './preds/'

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
  tst <- zero_check | sign_check
  return(tst)
}

#### 2. Predict ####
# obtain average percentage changes per EUNIS-ESy habitat lev.2 and 95%CI
res <- list() # results of summary stats
ab01 <- list() # results measuring percentages of plots with CMEIV change >= 0.1 and <= -0.1 
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
    mutate(eiv_abs.change = eiv_pred_max_yr - eiv_pred_min_yr) #%>%
    #select(-contains('_pred_'))
  
  # percentages of plots with CMEIV change <= -0.1 and >= 0.1
  ab01[[ind.name]] <- pred_focalchange %>%
    mutate(
      perc.below_e_0.1 = cut(eiv_abs.change, breaks = c(-Inf,-0.1,Inf), labels = c(1,0), include.lowest = F),
      perc.above_e_0.1 = cut(eiv_abs.change, breaks = c(-Inf,0.1,Inf), labels = c(0,1),  include.lowest = T),
    ) %>%
    gather('k','v',contains('0.1')) %>%
    group_by(ESy2, k,v) %>%
    summarise(
      n=n()
    ) %>%
    group_by(ESy2, k) %>%
    mutate(
      sum=sum(n)
    ) %>%
    filter(v==1) %>%
    select(-v) %>%
    mutate(
      percs = round((n/sum)*100, 2)
    ) %>%
    select(-n,-sum) %>%
    spread(k,percs) 
  
  # summarize stats
  res.i <- pred_focalchange %>%
    group_by(ESy2) %>%
    summarise(
      mean.min = mean(eiv_pred_min_yr),
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
  
  print(Sys.time()-st)
}

#### 3. Check homogenisation trends ####

homo_plts <- list()
for (ind.name in ind.names$eiv_name_raw) {
  homo_plts[[ind.name]] <- res[[ind.name]] %>%
    mutate(habitat=unname(sapply(.$ESy2, getESy1)), .before = ESy2) %>%
    ggplot(aes(mean.min, mean)) +
    geom_hline(yintercept = 0, lty=2) +
    geom_smooth(method = 'lm', color='black') +
    stat_cor(aes(label = paste( ..p.label.., sep = "~`,`~")),
             label.x = quantile(range(res[[ind.name]]$mean.min), .6),
             label.y = quantile(range(res[[ind.name]]$mean), .9),
             r.accuracy=0.01, p.accuracy = 0.001) +
    stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),
             label.x = quantile(range(res[[ind.name]]$mean.min), .6),
             label.y = quantile(range(res[[ind.name]]$mean), 1),
             r.accuracy=0.01, p.accuracy = 0.001) +
    geom_label(aes(label=ESy2, fill=habitat), size=2.5, alpha=.5, label.size=NA) +
    theme_bw() +
    labs(
      x = bquote(CM[EIV]~.(ind.names$eiv_name[ind.names$eiv_name_raw %in% ind.name])~" - baseline"),
      y = bquote(symbol('D')*CM[EIV]~.(ind.names$eiv_name[ind.names$eiv_name_raw %in% ind.name]))
    ) +
    theme(legend.position = 'none')
}

ggsave(filename=paste0(pth2preds, 'rtm.trends.raw.svg'),
       cowplot::plot_grid(homo_plts$EIV_L,
                          homo_plts$EIV_T,
                          homo_plts$EIV_M,
                          homo_plts$EIV_N,
                          homo_plts$EIV_R,
                          ncol = 2),
       height = 8, width = 6)
ggsave(filename=paste0(pth2preds, 'rtm.trends.raw.pdf'),
       cowplot::plot_grid(homo_plts$EIV_L,
                          homo_plts$EIV_T,
                          homo_plts$EIV_M,
                          homo_plts$EIV_N,
                          homo_plts$EIV_R,
                          ncol = 2),
       height = 8, width = 6)

#### 4. Format final table and export ####

res_spread <- res %>%
  map(\(x) {
    x %>%
      select(ESy2, EIV, mean) %>%
      spread(EIV, mean) %>%
      mutate(`Habitat type` = unname(sapply(.$ESy2, getESy1)), .before =
               ESy2) %>%
      arrange(`Habitat type`)
  })

res_n <- bind_rows(res) %>%
  group_by(ESy2) %>% summarise(n_median = median(n)) %>% ungroup() %>%
  mutate(`Habitat type` = unname(sapply(.$ESy2, getESy1)), .before = ESy2) %>%
  arrange(`Habitat type`)

res_fin <- res_n %>%
  left_join(res_spread$EIV_L, by = join_by(`Habitat type`, ESy2)) %>%
  left_join(res_spread$EIV_T, by = join_by(`Habitat type`, ESy2)) %>%
  left_join(res_spread$EIV_M, by = join_by(`Habitat type`, ESy2)) %>%
  left_join(res_spread$EIV_N, by = join_by(`Habitat type`, ESy2)) %>%
  left_join(res_spread$EIV_R, by = join_by(`Habitat type`, ESy2)) %>%
  left_join(hab_plot.size, by = join_by(ESy2)) %>%
  select(1, 2, plot_size_median, everything())

pretty_years <- dat.year.ESy2 %>%
  mutate(Years = paste0(min_yr, '-', max_yr)) %>%
  select(ESy2, Years)

res_fin <- res_fin %>%
  left_join(pretty_years, 'ESy2') %>%
  select(1:2, Years, everything())

# load ESy2 names table
hab.names <-
  read_delim('./data/EUNIS_ESy2_habitat.names.txt', col_names = F, show_col_types = F) %>%
  setNames(c('ESy2', 'Habitat name'))

res_fin <- res_fin %>%
  left_join(hab.names) %>%
  select(1:3, `Habitat name`, everything())

res_fin %>%
  select(Light:Reaction) %>%
  gather('k', 'v') %>%
  ggplot(aes(x = v)) +
  geom_histogram()

# export table with sign chars (without background color shades)
ft_dat <- res_fin %>%
  select(-ESy2) %>%
  rename(`no. plots` = n_median) %>%
  mutate(plot_size_median = round(plot_size_median, 0)) %>%
  mutate(
    Light = ifelse(Light > 0, paste0('+', Light), as.character(Light)) %>% str_pad(
      width = 5,
      side = c('right'),
      pad = '0'
    ),
    Temperature = ifelse(
      Temperature > 0,
      paste0('+', Temperature),
      as.character(Temperature)
    ) %>% str_pad(
      width = 5,
      side = c('right'),
      pad = '0'
    ),
    Moisture = ifelse(Moisture > 0, paste0('+', Moisture), as.character(Moisture)) %>% str_pad(
      width = 5,
      side = c('right'),
      pad = '0'
    ),
    Nitrogen = ifelse(Nitrogen > 0, paste0('+', Nitrogen), as.character(Nitrogen)) %>% str_pad(
      width = 5,
      side = c('right'),
      pad = '0'
    ),
    Reaction = ifelse(Reaction > 0, paste0('+', Reaction), as.character(Reaction)) %>% str_pad(
      width = 5,
      side = c('right'),
      pad = '0'
    )
  ) %>%
  rename(`Plot size` = plot_size_median)
ft_dat[ft_dat == '00000'] <- '0.00'

flextable(ft_dat) %>%
  save_as_docx(path = paste0(pth2preds, 'habitat_ESy2_means.docx'))

# export table with background color shades
ft_dat <- res_fin %>%
  select(-ESy2) %>%
  rename(`no. plots` = n_median) %>%
  mutate(plot_size_median = round(plot_size_median, 0)) %>%
  rename(`Plot size` = plot_size_median)

ft = ft_dat %>%
  flextable()

for (i in ind.names$eiv_name) {
  colpa <- c('#8395ce', '#d0e2de', '#ffffff', '#fac863', '#d15e6b')
  colorer <- col_bin(
    palette = colpa,
    bins = c(-Inf,-0.2,-0.1, 0.1, 0.2, Inf),
    domain = c(min(res_fin[[i]], na.rm = TRUE), max(res_fin[[i]], na.rm = TRUE)),
    # Provide the actual domain
    right = TRUE,
    reverse = F
  )
  ft <- ft %>%
    bg(j = i, bg = colorer, part = 'body')
}

ft %>%
  save_as_docx(path = paste0(pth2preds, 'habitat_ESy2_means_BGcol.docx'))

# Inspects percentages of plots with CMEIV change >= 0.1 and <= -0.1 
ab01_res <- ab01 %>%
  bind_rows(.id = 'eiv_name_raw') %>%
  left_join(ind.names) %>%
  left_join(hab.names) %>%
  left_join(res_n) %>%
  rename(`no. plots` = n_median) %>%
  mutate(`no. plots` = round(`no. plots`)) %>%
  ungroup() 
a01 <- ab01_res %>%
  select(`Habitat type`, `Habitat name`, `no. plots`, eiv_name, perc.above_e_0.1) %>%
  spread(4, 5)
b01 <- ab01_res %>%
  select(`Habitat type`, `Habitat name`, `no. plots`, eiv_name, perc.below_e_0.1) %>%
  spread(4, 5)
a01 %>%
  flextable() %>%
  save_as_docx(path = paste0(pth2preds, 'habitat_ESy2_perc.above_e_0.1.docx'))
b01 %>%
  flextable() %>%
  save_as_docx(path = paste0(pth2preds, 'habitat_ESy2_perc.below_e_0.1.docx'))

# quit
quit(save = 'no')