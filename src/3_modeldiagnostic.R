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
  eiv_name_raw = c('EIV_M', 'EIV_N', 'EIV_T', 'EIV_L', 'EIV_R'),
  eiv_name = c(
    'Moisture',
    'Nutrients',
    'Temperature',
    'Light',
    'Soil reaction'
  )
)

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

# import last fit models
m <- ind.names$eiv_name_raw %>%
  map(\(ind.name) {
    list.files(
      './models',
      pattern = paste0('RF.last_fit_', ind.name),
      full.names = T
    ) %>%
      read_rds()
  }) %>%
  setNames(ind.names$eiv_name_raw)

# import CV results
cv_res <- ind.names$eiv_name_raw %>%
  map(\(ind.name) {
    list.files(
      './models',
      pattern = paste0('RF.cv_res_', ind.name),
      full.names = T
    ) %>%
      read_rds()
  }) %>%
  setNames(ind.names$eiv_name_raw)

# import tuning results
tune_res <- ind.names$eiv_name_raw %>%
  map(\(ind.name) {
    list.files(
      './models',
      pattern = paste0('RF.tune_res_', ind.name),
      full.names = T
    ) %>%
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

#### 4. Last fit evaluation metrics ####\
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
    (
      extract_workflow(x) %>%
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
  varimpres$EIV_M,
  varimpres$EIV_N,
  varimpres$EIV_T,
  varimpres$EIV_L,
  varimpres$EIV_R
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
  geom_abline(
    lty = 2,
    color = 'red',
    lwd = .5,
    alpha = .8
  ) +
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
# get EU sf data
library(rnaturalearth)
regions_name <- c('Albania', 'Austria', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria',
                  'Corsica', 'Crete', 'Croatia', 'Czechia', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany',
                  'Greece', 'Hungary', 'Ireland', 'Italy', 'Kosovo', 'Latvia', 'Liechtenstein',
                  'Lithuania', 'Luxembourg', 'Malta', 'Moldova', 'Montenegro', 'Netherlands', 'North Macedonia',
                  'Norway', 'Poland', 'Portugal', 'Romania', 'Sardinia', 'Serbia', 'Sicily',
                  'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine', 'United Kingdom')
bbox_coords <- c(xmin = -1123055, ymin = 3923814, xmax = 2796649, ymax = 8007282)
EU <- ne_countries(scale = 'large', returnclass = 'sf') %>%
  filter(name %in% regions.name) %>%
  st_transform(crs = 25832) %>%
  st_crop(bbox_coords) %>% 
  select(geometry)

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
dip <- list()
for (i in c('EIV_M', 'EIV_N', 'EIV_T', 'EIV_L', 'EIV_R')) {
  dip[[i]] <- paste0('./preds/',
                     i,
                     '.preds.rf.csv.gz') %>%
    read_csv(show_col_types = F) %>%
    select(plot_id, habitat, contains('1960.2020'), lm.slope_estimate) %>%
    mutate(eiv_name_raw = i)
}

eiv_abs.change_data <- list()
for (i in c('EIV_M', 'EIV_N', 'EIV_T', 'EIV_L', 'EIV_R')) {
  eiv_abs.change_data[[i]] <- dip[[i]] %>%
    select(plot_id, habitat, eiv_abs.change_1960.2020) %>%
    setNames(c('plot_id', 'habitat', i))
}
eiv_abs.change_data <- reduce(eiv_abs.change_data, left_join, by = c('plot_id', 'habitat'))
head(eiv_abs.change_data)

getsamplesize <- function(vec_a, vec_b) {
  nrow(drop_na(data.frame(vec_a, vec_b)))
}
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}
plot.inc.cor <- function(mat,
                         round.label.digit = 2,
                         cor.method = 'pearson',
                         size.label = 4) {
  corm <- cor(mat, method = cor.method, use = "pairwise.complete.obs")
  
  pcorm <- cor(mat, method = cor.method, use = "pairwise.complete.obs") %>%
    get_lower_tri %>%
    reshape2::melt()
  pcorm$value = ifelse(pcorm$Var1 == pcorm$Var2, NA, pcorm$value)
  psams <- corrr::colpair_map(mat, getsamplesize) %>%
    as.data.frame %>%
    tibble::column_to_rownames('term') %>%
    as.matrix %>%
    get_upper_tri %>%
    reshape2::melt() %>%
    setNames(names(pcorm))
  
  psams$value = ifelse(is.na(psams$value),
                       NA,
                       prettyNum(psams$value, big.mark = ',', scientific = F))
  
  pcorm$value = round(pcorm$value, round.label.digit)
  # Format the numeric values with two decimal places
  pcorm$formatted_value <- sprintf(paste0('%.', round.label.digit, 'f'), pcorm$value)
  pcorm$formatted_value <- ifelse(pcorm$formatted_value == 'NA', NA, pcorm$formatted_value)
  
  p <- ggplot() +
    geom_tile(data = pcorm,
              aes(x = Var1, y = Var2, fill = value),
              color = 'grey80') +
    geom_text(
      data = psams,
      aes(x = Var1, y = Var2, label = value),
      color = "black",
      size = size.label
    ) +
    geom_text(
      data = pcorm,
      aes(x = Var1, y = Var2, label = formatted_value),
      color = "black",
      size = size.label
    ) +
    scale_fill_gradient2(
      low = "#4A6FE3",
      mid = "white",
      high = "#D33F6A",
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab",
      na.value = 'white',
      name = paste0(firstup(cor.method), '\ncorrelation')
    ) +
    theme_classic() +
    guides(alpha = 'none') +
    theme(
      axis.text.x = element_text(
        size = 12,
        angle = 45,
        vjust = 1,
        hjust = 1,
        color = 'black'
      ),
      axis.text.y = element_text(size = 12, color = 'black'),
      axis.title = element_blank(),
      axis.line = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white")
    ) +
    coord_fixed() +
    geom_abline(slope = 1,
                intercept = 0,
                color = 'grey80')
  return(p)
}

plot_commlevel <- list()
for (h in c('Forest', 'Grassland', 'Scrub', 'Wetland')) {
  plot_commlevel[[h]] <- eiv_abs.change_data %>%
    filter(habitat == h) %>%
    select(contains('EIV')) %>%
    setNames(c('Moist.', 'Nutr.', 'Temp.', 'Light', 'Reac.')) %>%
    plot.inc.cor() +
    ggtitle(h) +
    theme(title = element_text(face = 2))
}

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