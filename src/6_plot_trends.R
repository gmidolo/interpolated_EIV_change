################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 09.07.2025
################################################################################

# Description: Get summary figure for EIV trends (1960-2020) across habitats

################################################################################


# 1. prepare data####

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# set indicator names
ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

# source help functions 
source('./src/0_helpfunctions.R')

# folder where to store figures and diagnostic output
pth2fig <- './fig/'

# load EVA and ReSurveyEurope data
set.seed(123)
dat.initial <-
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope( # Load ReSurveyEU (static, using one random point in the survey)
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
    )[['traintest_data']]  
  )

# load predictions for each EIV variable (needed for histogram plotting)
dat <- ind.names$eiv_name_raw %>%
  map(\(ind.name){
    # get plot id to discard
    th <- dat.initial[,c('plot_id', paste0('n.', ind.name))] %>%
      setNames(c('plot_id','th')) %>%
      filter(th < 0.8)
    # get predictions
    y <- paste0('./preds/', ind.name, '.preds.rf.csv.gz') %>%
      read_csv(show_col_types = F) %>%
      filter(year >= 1960 & year <= 2020) %>%
      rename(eiv_name_raw = eiv_name) %>%
      left_join(ind.names, 'eiv_name_raw') %>%
      select(-eiv_name_raw) %>%
      select(eiv_name, everything()) %>%
      arrange(plot_id) %>%
      anti_join(th, 'plot_id')
    return(y)
  }) %>%
  setNames(ind.names$eiv_name_raw)

# get no. plots available per habitat and eiv variable
plots_count <- dat %>% 
  bind_rows(.id='eiv_name_raw') %>%
  group_by(eiv_name_raw, habitat) %>%
  summarise(n=n(), .groups = 'drop')

# 2. Plot trends over time per habitat####

# load summary stats of predictions from 1960 to 2020
dat_alltrees <- ind.names$eiv_name_raw %>%
  map(\(ind.name){
    paste0('./preds/pdp/', ind.name, '_habitat_pdp_trends_alltrees.csv') %>%
      read_csv(show_col_types = F) %>%
      mutate(eiv_name_raw=ind.name) %>%
      left_join(ind.names, by = 'eiv_name_raw') %>%
      select(contains('eiv_name'), everything())
  }) %>%
  setNames(ind.names$eiv_name_raw)
  
# calculate slope across studied period
slope_trend <- dat_alltrees %>%
  map(\(x){
    x %>% 
      group_by(habitat) %>%
      do(broom::tidy(lm(mean ~ year , data = .))) %>%
      filter(term == 'year')
  }) %>% 
  bind_rows(.id='eiv_name_raw') %>%
  mutate(slope_decades = estimate*10)
slope_trend$slope_decades %>% hist()
slope_trend$slope_decades %>% range()

# generate color palettes over linear slope changes
vals = unique(round(seq(-0.06, 0.06, 0.001), 2))
pals = hcl.colors(length(vals), palette = 'RdYlBu', rev=T)

## generate and export legend
# pdf(paste0(pth2fig,'legenda_average.trends.pdf'), height = 5, width = 3.5)
# legend_image <- as.raster(matrix(rev(colorspace::adjust_transparency(pals, 0.95)), ncol=1))
# plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'EIV gained/loss\nper decade')
# text(x=1.5, y = seq(0,1,l=5), labels = quantile(as.numeric(vals), seq(0,1,.25)))
# rasterImage(legend_image, 0, 0, 1,1)
# dev.off()

# generate various baseline EIV or average EIV per year 
dat_cntrd <- dat_alltrees %>%
  map(\(x){

    dat_baseline <- x %>% filter(year==1960) 
    names(dat_baseline)[5:length(dat_baseline)] <- paste0(names(dat_baseline)[5:length(dat_baseline)],'_1960')
    
    dat_cntrd_x <- x %>%
      left_join(
        dat_baseline %>% select(-year), by='habitat'
      ) %>%
      group_by(habitat) %>%
      mutate(
        m = mean-mean_1960,
        lpi = m + (quantile0.025 - median), # lower prediction interval - quantile
        upi = m + (quantile0.975 - median), # upper prediction interval - quantile
        lci = m + (lower_ci_95 - mean), # lower confidence interval of the mean
        uci = m + (upper_ci_95 - mean) # upper confidence interval of the mean
      ) %>%
      select(habitat, year, mean_1960, m, lci, uci, lpi, upi)
    
    return(dat_cntrd_x)
  })

# Plot figure (quick figure)
dat_cntrd %>%
  bind_rows(.id='eiv_name_raw') %>%
  ggplot(aes(year, m))+
  geom_ribbon(aes(ymin = lpi, ymax = upi), fill = 'grey70')+
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = 'grey30')+
  geom_line()+
  facet_wrap(eiv_name_raw~habitat, scales='free_y')

# ranges in y-limits EIV (prediction interval)
range(
  dat_cntrd %>% bind_rows(.id='eiv_name_raw') %>% pull(lpi),
  dat_cntrd %>% bind_rows(.id='eiv_name_raw') %>% pull(upi)
)

# define limits on y axis
y_breaks <- seq(-0.6, 0.6, 0.3)
y_limits_centrd <- range(y_breaks)

# plot trends for each EIV-habitat combination
p_means.1960 <- list()

hab <- c('Forest','Grassland','Scrub','Wetland')

for(i in names(dat_cntrd)) {
  
  for (h in hab) {
  
  # define background color based on the slope estimate
  slope_hi <- slope_trend %>% 
    filter(habitat==h & eiv_name_raw == i) %>% 
    pull(slope_decades)
  bg_col <- pals[which(vals %in% round(slope_hi, 2))]
  
  # get plot count to display
  nplots <- plots_count %>%
    filter(habitat==h & eiv_name_raw == i) %>%
    prettyNum(big.mark=",",scientific=FALSE)

  # define y axis lables 
  y_left <- y_breaks[1:floor(length(y_breaks)/2)]
  y_right <- y_breaks[(ceiling(length(y_breaks)/2)+1):length(y_breaks)]
  y_right <- paste0('+', y_right)
  y_mid <- dat_cntrd[[i]] %>% 
            filter(habitat==h) %>% 
            pull(mean_1960) %>% 
            unique() %>% 
            round(1)
  y_labels <- c(y_left, y_mid, y_right)

  # get pretty EIV variable name
  pretty_EIV_name <- ind.names[which(ind.names$eiv_name_raw==i),]$eiv_name
  
  # plot relative changes over 1960 as baseline
  p_means.1960[[h]][[pretty_EIV_name]] <-
    dat_cntrd[[i]] %>%
    filter(habitat==h) %>%
      ggplot(aes(year, m)) +
      # color background
      annotate('rect', xmin=1960, xmax=2020, ymin=-Inf, ymax=Inf, alpha=.8, fill=bg_col) +
      # plot preidction interval
      geom_ribbon(aes(ymin=lpi, ymax=upi), fill='grey80', alpha=.7) +
      # plot CI to the mean
      geom_ribbon(aes(ymin=lci, ymax=uci), fill='grey20', alpha=.7) +
      # plot zero change line
      geom_segment(x = 1960, y = 0, xend = 2020, yend = 0, lty=3, color='grey30', linewidth=.3) +
      # details
      labs(y=paste0('Mean EIV ', pretty_EIV_name), x='Year') +
      theme_bw()+
      theme(axis.title = element_blank(),
            panel.grid = element_blank())+
      # set appropriate scale 
      scale_y_continuous(
        limits = y_limits_centrd,
        labels = y_labels,
        breaks = y_breaks
    ) +
    #annotate no. plots
    annotate('text', x = 1990, y = Inf, label = paste0('n = ', nplots['n']), vjust = 1.5, size= 2.75)
 }
}

# 3. Plot histogram-like column bars####

p.hists <- list()

# define breaks & labels
br <- c(-Inf,-1, -0.5, -0.1, 0.1, 0.5, 1, Inf)
lb <- c('< 1','-1 - -0.5', '-0.5 - -0.1', '-0.1 - 0.1', '0.1 - 0.5', '0.5 - 1', '> 1')

# define color palette
colpa <- data.frame(cols=hcl.colors(length(lb), palette = 'RdYlBu', rev=T), mean_cat=lb)

# categorize data into groups
dat_cat <- dat %>% 
  map(\(x){
    x %>%
     mutate(mean_cat = cut(
     eiv_abs.change_1960.2020, breaks = br, labels = lb)) %>% 
     group_by(mean_cat) %>%
     summarise(n = n())
  }) 

dat_mean.cat <- dat %>% 
  map(\(x){
   y=as.character(round(mean(x$eiv_abs.change_1960.2020), 2))
   if(substr(y, start = 1, stop = 1)!='-'){
     y=paste0('+',y)
   }
   return(y)
  }) 

dat_mean.val <- dat %>% 
  map(\(x){
   round(mean(x$eiv_abs.change_1960.2020), 4)
  }) 

# plot the histograms
for(i in names(dat)) { 
  di <- dat_cat[[i]] 
  di$count <- 1:nrow(di)
  # get pretty EIV variable name
  pretty_EIV_name <- ind.names[which(ind.names$eiv_name_raw==i),]$eiv_name
  norm.mean <- weighted.mean(di$count, di$n)
  max_count <- di %>% filter(n == max(n)) %>% pull(n)  
  txt_y_pos <- max_count-max_count*0.05
  p.hists[[pretty_EIV_name]] <- di %>%
      ggplot(aes(x = count, y = n, fill =  mean_cat)) +
      geom_col(col=NA, position = 'dodge', width=1) +
      geom_vline(xintercept = mean(di$count), lty=2, col='grey50') +
      geom_vline(xintercept = norm.mean, lty=1) +
      annotate('text', x = 6, y = txt_y_pos, label=dat_mean.cat[[i]]) +
      scale_fill_manual(values=colpa$cols) +
      theme_bw() +
      ggtitle(pretty_EIV_name)+
      scale_y_continuous(labels = scales::unit_format(unit = '', scale = 1e-3))+
      theme(legend.position = 'none', 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.title = element_blank(),
            panel.border = element_blank(), 
            axis.line.y = element_line(colour = 'black'),
            axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# combine plots into a grid
cp <- plot_grid(
  # Row 1: Histograms
  p.hists$Light,
  p.hists$Temperature,
  p.hists$Moisture,
  p.hists$Nitrogen,   
  p.hists$Reaction,
  
  # Row 2: Forest means
  p_means.1960$Forest$Light,
  p_means.1960$Forest$Temperature,
  p_means.1960$Forest$Moisture,
  p_means.1960$Forest$Nitrogen,
  p_means.1960$Forest$Reaction,
  
  # Row 3: Grassland means
  p_means.1960$Grassland$Light,
  p_means.1960$Grassland$Temperature,
  p_means.1960$Grassland$Moisture,
  p_means.1960$Grassland$Nitrogen,
  p_means.1960$Grassland$Reaction,
  
  # Row 4: Scrub means
  p_means.1960$Scrub$Light,
  p_means.1960$Scrub$Temperature,
  p_means.1960$Scrub$Moisture,
  p_means.1960$Scrub$Nitrogen,
  p_means.1960$Scrub$Reaction,
  
  # Row 5: Wetland means
  p_means.1960$Wetland$Light,
  p_means.1960$Wetland$Temperature,
  p_means.1960$Wetland$Moisture,
  p_means.1960$Wetland$Nitrogen,
  p_means.1960$Wetland$Reaction,
  
  cols = 5,
  rows = 5
)

cp

# export figure
ggsave(paste0(pth2fig, 'minitrends.changeto1960.jpg'), cp, width = 8, height = 7.5, dpi = 600)
# ggsave(paste0(pth2fig, 'minitrends.changeto1960.pdf'), cp, width = 8, height = 7.5)
