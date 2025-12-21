################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Diagnostic / evaluation of Random Forests models

################################################################################

#### 1. Data preparation / import ####

# define pretty names for EIV indicators (to be used in plots)
ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

# source help functions 
source('./src/0_helpfunctions.R')

# folder with model results
pth2export <- './models/'

# folder where to store figures and diagnostic output
pth2fig <- './fig/'

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(flextable)
  library(vip)
  library(sf)
  library(terra)
})

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

# import last fit models
m <- ind.names$eiv_name_raw %>%
  map(\(ind.name) {
    paste0(pth2export, 'RF.last_fit_', ind.name, '.rds') %>%
      read_rds()
  }) %>%
  setNames(ind.names$eiv_name_raw)

# import CV results
cv_res <- ind.names$eiv_name_raw %>%
  map(\(ind.name) {
    paste0(pth2export, 'RF.cv_metrics_', ind.name, '.csv') %>%
      read_csv(show_col_types = F)
  }) %>%
  setNames(ind.names$eiv_name_raw)

# import tuning results
tune_res <- ind.names$eiv_name_raw %>%
  map(\(ind.name) {
    paste0(pth2export, 'RF.tune_res_', ind.name, '.rds') %>%
      read_rds()
  }) %>%
  setNames(ind.names$eiv_name_raw)

#### 2. Plot tuning results ####
# prepare data
tune_res_tidy <- tune_res %>%
  map(\(x) {
    (ggplot_build(autoplot(x)))$plot$data
  }) %>%
  bind_rows(.id = 'eiv_name_raw') %>%
  left_join(ind.names, 'eiv_name_raw')
tune_res_tidy$`Minimal Node Size` <- as.factor(tune_res_tidy$`Minimal Node Size`)
tune_res_tidy$.metric <- ifelse(tune_res_tidy$.metric == 'rmse', 'RMSE', 'R-squared')
tune_res_tidy$eiv_name <- factor(tune_res_tidy$eiv_name, ind.names$eiv_name)
best_res_tidy <- tune_res_tidy %>% filter(.metric == 'RMSE') %>% group_by(eiv_name) %>% filter(mean == min(mean))
best_res_tidy <- bind_rows(
  best_res_tidy,
  tune_res_tidy %>% filter(.metric == 'R-squared') %>% semi_join(
    best_res_tidy,
    by = c('Minimal Node Size', 'value', 'eiv_name_raw')
  )
)
best_res_tidy$eiv_name <- factor(best_res_tidy$eiv_name, ind.names$eiv_name)

# plot
tune_res_plot <- ggplot(tune_res_tidy, aes(x = value, y = mean, col = `Minimal Node Size`)) +
  geom_point() +
  geom_line() +
  facet_wrap(.metric ~ eiv_name, scales = 'free_y', nrow = 2) +
  scale_color_manual(values = hcl.colors(5, 'Blues', rev = T)[-1]) +
  geom_point(data = best_res_tidy, aes(x = value, y = mean), col = 'red') +
  theme_bw() +
  theme(axis.title.y = element_blank(), legend.position = 'top') +
  labs(x = 'mtry') +
  scale_y_continuous()
tune_res_plot

# export
ggsave(
  paste0(pth2fig, 'RF.tune_res_plot.png'),
  tune_res_plot,
  width = 10,
  height = 4,
  dpi = 600
)


#### 3. 10-fold CV results ####
# prepare data
cv_res_tidy <- cv_res %>%
  map(collect_metrics) %>%
  bind_rows(.id = 'eiv_name_raw') %>%
  left_join(ind.names, 'eiv_name_raw') %>%
  select(eiv_name, .metric, mean, std_err)
cv_res_tidy$.metric <- ifelse(cv_res_tidy$.metric == 'rmse', 'RMSE', 'R-squared')
cv_res_tidy

# aggregate results & export
cv_res_tidy %>%
  mutate(std_err = round(std_err, 6), mean = round(mean, 3)) %>%
  setNames(c('EIV', 'metric', 'mean', 'std. err')) %>%
  flextable() %>%
  merge_v('EIV') %>%
  autofit() %>%
  save_as_docx(path = paste0(pth2fig, 'RF.cv_res.docx'))

#### 4. Last fit evaluation metrics ####
# eval each model over testing data
m %>% 
  map(collect_metrics) 

#### 5. Display variable importance #### 
# set pretty names to display
pretty_vars <- data.frame(
  Variable = c('habitat', 'y', 'x', 'elev', 'n', 'year', 'plot_size'),
  Variable_name = c(
    'Habitat',
    'Northing',
    'Easting',
    'Elevation',
    'No. Species',
    'Time',
    'Plot Size'
  )
)

# retrieve var. importance via vip::vip()
var_imp <- m %>%
  map(\(x) {
    (extract_workflow(x) %>%
        extract_fit_parsnip() %>%
        vip(geom = 'col') %>%
        ggplot_build()
    )$plot$data %>%
      left_join(pretty_vars, 'Variable')
  }) %>%
  bind_rows(.id = 'eiv_name_raw') %>%
  left_join(ind.names, 'eiv_name_raw') %>%
  mutate(eiv_name = factor(eiv_name, ind.names$eiv_name))

# plot
varimpres <- list()
for (i in ind.names$eiv_name_raw) {
  di <- var_imp %>% filter(eiv_name_raw == i) %>% arrange(Importance)
  di$Variable_name <- factor(di$Variable_name, di$Variable_name)
  varimpres[[i]] <- ggplot(di, aes(x = Variable_name, y = Importance)) +
    geom_col() +
    coord_flip() +
    ggtitle(unique(di$eiv_name)) +
    theme_bw() +
    theme(axis.title.y = element_blank())
}

# aggregate plots
varimpres_plot <- cowplot::plot_grid(
  varimpres$EIV_L,varimpres$EIV_T,varimpres$EIV_M,varimpres$EIV_N,varimpres$EIV_R
)

# export
ggsave(
  paste0(pth2fig, 'RF.varimp.jpg'),
  varimpres_plot,
  width = 10,
  height = 6,
  dpi = 600
)

#### 6. Plot last fit  #### 
# collect predictions over testing set
set.seed(123)
pred_test <- ind.names$eiv_name_raw %>% map(\(x) {
  cat(x, ' |> ') # show progress
  d <- dat_split
  names(d)[which(names(d) == paste0('n.', x))] <- 'treshold'
  names(d)[which(names(d) == paste0('cm.', x))] <- 'eiv'
  treshold.of.EIVE.species = 0.8 
  d <- d %>%
    filter(treshold >= treshold.of.EIVE.species) %>%
    select(plot_id, eiv, x, y, elev, year, habitat, n, plot_size) %>%
    initial_split(prop = 4/5, strata = eiv)
  tstd <- testing(d)
  agu_res <- augment(extract_workflow(m[[x]]), tstd)
  return(agu_res)
}) %>%
  setNames(ind.names$eiv_name_raw) %>%
  bind_rows(.id = 'eiv_name_raw') %>%
  left_join(ind.names, 'eiv_name_raw') %>%
  mutate(eiv_name = factor(eiv_name, ind.names$eiv_name))

# define set of metrics to evaluate the model
eval_metrics <- metric_set(rmse, rsq, rpd)

# eval over obs. vs. pred. data 
pred_test_res <- pred_test %>%
  group_split(eiv_name) %>%
  map(\(x) {
    x %>%
      eval_metrics(eiv, .pred)
  }) %>%
  setNames(ind.names$eiv_name_raw) %>%
  bind_rows(.id = 'eiv_name_raw') %>%
  left_join(ind.names, 'eiv_name_raw') %>%
  mutate(eiv_name = factor(eiv_name, ind.names$eiv_name)) %>%
  select(eiv_name, .metric, .estimate) %>%
  spread(2, 3) %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(txt = paste0('RMSE = ', rmse, '\nRSQ = ', rsq, '\nRPD = ', rpd)) %>%
  left_join(
    pred_test %>%
      group_by(eiv_name) %>%
      summarise(.pred = quantile(range(.pred), 0.9), eiv =   quantile(range(eiv)  , 0.18)),
    'eiv_name'
  )

# plot
pred_test_plot <- ggplot(pred_test, aes(eiv, .pred)) +
  geom_hex(bins = 45) +
  scale_fill_gradient(
    low = 'lightblue',
    high = 'midnightblue',
    trans = 'log10'
  ) +
  labs(x = expression(paste('Observed ', CM[EIV], ' (test set)')),
       y = expression(paste('Predicted ', CM[EIV])),
       fill = 'No. plots') +
  geom_abline(lty = 2, color = 'red', lwd = .5, alpha = .8) +
  theme_bw() + theme(legend.position = 'bottom') +
  facet_wrap( ~ eiv_name, scales = 'free') +
  geom_text(
    data = pred_test_res,
    aes(label = txt),
    fill = NA,
    alpha = .8,
    size = 3
  )

# export
ggsave(
  paste0(pth2fig, 'RF.pred_test_plot.png'),
  pred_test_plot,
  width = 9 * 0.85,
  height = 7.5 * 0.85,
  dpi = 600
)


#### 7. Geographic patterns in model residuals #### 
# get EU-countires sf shapes
EU <- read_rds('./data/EU_shape_map.rds') %>% 
  st_buffer(1000) %>% 
  st_simplify(dTolerance = 4000)

# calculate residuals (= obs - pred)
pred_test$resid <- (pred_test$eiv - pred_test$.pred) # model residuals
hist(pred_test$resid)

# define categorical levels for plotting 
categorize_rasters = TRUE
if (categorize_rasters) {
  br = c(-Inf, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, Inf)
  lb = c('< -0.2',
         '-0.2 - -0.1',
         '-0.1 - -0.05',
         '-0.05 - 0.05',
         '0.05 - 0.1',
         '0.1 - 0.2',
         '> 0.2')
  df.cols = data.frame(lb, cols = hcl.colors(length(lb), 'RdYlBu', rev = T))
}

# aggregate (average) residuals value into 50 km rasters
rasters <- list()
contas <- list()
for (i in ind.names$eiv_name_raw) {
  for (h in c('Forest', 'Grassland', 'Scrub', 'Wetland')) {
    pdh <- pred_test %>% filter(habitat == h & eiv_name_raw == i)
    contas[[i]][[h]] <- pdh %>%
      summarise(n = paste0('n = ', prettyNum(
        n(), big.mark = ',', scientific = F
      )))
    
    min5 <- \(x){sum(x) >= 5} # at least 5 plots available?
    r0 <- rast(res = 50 * 1000, # base (empty) raster (res in meters)
               extent = ext(EU),
               crs = crs(EU))
    r.min5 <- rasterize(pdh %>% select(x, y) %>% as.matrix(),
                        r0,
                        values = 1,
                        fun = min5)
    r.min5 <- clamp(r.min5, lower = 1, value = FALSE)
    r <- rasterize(pdh %>% select(x, y) %>% as.matrix(),
                   r0,
                   values = pdh$resid,
                   fun = mean)
    r <- mask(r, r.min5)
    newvals.cont = as.vector(values(r))
    newvals.cat = cut(newvals.cont, breaks = br, labels = lb)
    if (categorize_rasters) {
      values(r) <- newvals.cat
    } else {
      values(r) <- newvals.cont
    }
    rasters[[i]][[h]] <- r
  }
}

# bind raster data
rasters_data <- rasters %>%
  map(\(x){
    x %>% map(\(k) {
      as.data.frame(k, xy = T)
    }) %>% bind_rows(.id = 'habitat')
  }) %>% bind_rows(.id = 'eiv_name_raw') %>%
  left_join(ind.names, 'eiv_name_raw') %>%
  mutate(eiv_name = factor(eiv_name, ind.names$eiv_name))

# bind counts (no. plots)
contas_data <- contas %>%
  map(\(x){
    x %>% map(\(k) {
      as.data.frame(k)
    }) %>% bind_rows(.id = 'habitat')
  }) %>% bind_rows(.id = 'eiv_name_raw') %>%
  left_join(ind.names, 'eiv_name_raw') %>%
  mutate(eiv_name = factor(eiv_name, ind.names$eiv_name))

# plot maps
map_residuals <- ggplot() +
  geom_sf(data = EU, col = NA) +
  geom_raster(data = rasters_data, aes(x, y, fill = mean)) +
  theme(axis.title = element_blank()) +
  labs(fill = expression(paste(CM[EIV], ' residuals'))) +
  scale_fill_manual(values = df.cols$cols) +
  geom_text(
    data = contas_data,
    x = -263716.1,
    y = 7629361,
    aes(label = n),
    fontface = 1,
    size = 3.25
  ) +
  facet_grid(eiv_name ~ habitat) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.position = 'right',
    strip.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )

# export map
ggsave(
  paste0(pth2fig, 'RF.map_distribution_of_residuals.png'),
  map_residuals,
  width = 9,
  height = 13,
  dpi = 600
)

#### 8. Correlations across predictions of change####
# N.B. this part requires having collected predictions (using `./src/interpolation.R`)
set.seed(123)
dat.initial <-
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope( # Load ReSurveyEU (static, using one random point in the survey)
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
    )[['traintest_data']]  
  )

dat <- ind.names$eiv_name_raw %>%
  map(\(ind.name){
    # get plot id to discard
    th <- dat.initial[,c('plot_id', paste0('n.', ind.name))] %>%
      setNames(c('plot_id','th')) %>%
      filter(th<0.8)
    # get predictions
    y <- paste0('./preds/', ind.name, '.preds.rf.csv.gz') %>%
      read_csv(show_col_types = F) %>%
      filter(year >= 1960 & year <= 2020) %>%
      select(plot_id, habitat, eiv_abs.change_1960.2020) %>%
      arrange(plot_id) %>%
      anti_join(th, 'plot_id')
    return(y)
  })

# set new names
names(dat) <- ind.names$eiv_name
for (i in names(dat)) {
  dat[[i]] <- dat[[i]] %>% setNames(c('plot_id','habitat',i))
}

# reduce list of tibbles with a full join
dat <- reduce(dat, full_join, by=c('plot_id','habitat'))
nrow(dat)

# plot
plot_commlevel <- list()
for (h in c('Forest', 'Grassland', 'Scrub', 'Wetland')) {
  plot_commlevel[[h]] <- dat %>%
    filter(habitat == h) %>%
    select(-plot_id, -habitat) %>%
    setNames(c('Light', 'Tempe.', 'Moist.', 'Nitro.', 'React.')) %>%
    plot.inc.cor() +
    ggtitle(h) +
    theme(title = element_text(face = 2))
}

# combine plot
cp <- cowplot::plot_grid(
  plot_commlevel$Forest,
  plot_commlevel$Grassland,
  plot_commlevel$Scrub,
  plot_commlevel$Wetland
)
cp

ggsave(
  paste0(pth2fig, 'corrplot_CM.change_raw.pdf'),
  cp,
  width = 10,
  height = 8
)


## quit ##
quit(save = 'no')