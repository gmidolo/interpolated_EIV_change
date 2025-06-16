################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.05.2025
################################################################################

# Description: Get summary figure for EIV trends (1960-2020) across habitats

################################################################################

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

# folder where data with predictions (in .csv) were exported
pth2preds <- './preds/'

# folder where to store figures and diagnostic output
pth2fig <- './fig/'

# prepare data for interpolation
set.seed(123)
dat.initial <-
  bind_rows(
    read_csv('./data/EVA.csv.xz', show_col_types = F), # Load EVA data
    format_ReSurveyEurope( # Load ReSurveyEU (static, using one random point in the survey)
      training_strategy = 'random',
      path_resurvey_clean = './data/ReSurveyEU_clean.csv.xz'
    )[['traintest_data']]  
  )

# load predictions for each EIV variable
dat <- ind.names$eiv_name_raw %>%
  map(\(ind.name){
    # get plot id to discard
    th <- dat.initial[,c('plot_id', paste0('n.', ind.name))] %>%
      setNames(c('plot_id','th')) %>%
      filter(th<0.8)
    # get predictions
    y <- paste0(pth2preds, ind.name, '.preds.rf.csv.gz') %>%
      read_csv(show_col_types = F) %>%
      filter(year >= 1960 & year <= 2020) %>%
      rename(eiv_name_raw = eiv_name) %>%
      left_join(ind.names, 'eiv_name_raw') %>%
      select(-eiv_name_raw) %>%
      select(eiv_name, everything()) %>%
      arrange(plot_id) %>%
      anti_join(th, 'plot_id')
    return(y)
  })

# set new names
names(dat) <- ind.names$eiv_name

# get plot counts
dat %>% bind_rows() %>% select(plot_id) %>% unique() %>% nrow()

# get general trends everywhere
dat %>% map(\(x){mean(x$eiv_perc.change_1960.2020)})

# get general trends by ESy lev-1
dat %>% map(\(x){x %>% group_by(habitat) %>% summarise(mean.perc=mean(eiv_perc.change_1960.2020))})
  
### define periods
periods <- c('From 1960 to 1980', 'From 1980 to 2000', 'From 2000 to 2020', 'From 1960 to 2020')
period.dat <- data.frame(
  periods=periods,
  min=c(1960, 1980, 2000, 1960),
  max=c(1980, 2000, 2020, 2020)
)

# no. plots used for interpolation
dat %>% bind_rows() %>% select(plot_id, year) %>% filter(year >= 1960 & year <= 2020) %>% select(plot_id) %>% unique %>% nrow

# calc mean values
means_dats <- dat %>%
  map(\(x){
    x %>% 
    filter(year >= 1960 & year <= 2020) %>% # select only plots sampled within that period
    select(plot_id, eiv_name, habitat, contains('eiv_pred_')) %>%
    gather('year', 'eiv_pred', contains('eiv_pred_')) %>%
    mutate(year = as.numeric(gsub("\\D", "", year))) %>%
    group_by(habitat, year) %>%
    summarise(mean = mean(eiv_pred), sd = sd(eiv_pred), n=n())  %>%
    ungroup() 
  }) 

# calculate slope across studied period(s)
period_slope = list()
period = list()
for(k in names(means_dats)) {
  message('Starting trend estimation for: ', k)
  for (i in periods[4]) { # We will just focus on the entire period (from 1960 to 2020)
    st=Sys.time()
    cat('Period: ', i,' - Start')
    
    max_i = period.dat %>% filter(periods==i) %>% pull(max)
    min_i = period.dat %>% filter(periods==i) %>% pull(min)
    
    period[[k]][[i]] <- dat[[k]] %>%
      select(plot_id, habitat, contains('eiv_pred_')) %>%
      gather('year', 'eiv_pred', contains('eiv_pred_')) %>%
      mutate(year = as.numeric(gsub("\\D", "", year))) %>%
      filter(year >= min_i & year <= max_i) 
        
    # calc slope estimate
    period_slope[[k]][[i]] <- period[[k]][[i]]%>%
      group_by(habitat) %>%
      do(broom::tidy(lm(eiv_pred ~ year , data = .))) %>%
      filter(term == 'year')
    cat('Period: ', i,' - Done')
    print(Sys.time()-st)
    
  }
  period_slope[[k]] <- period_slope[[k]] %>% bind_rows(.id='period')
  period[[k]] <- period[[k]] %>% bind_rows(.id='period')

  period_slope[[k]]$slope_decades <- period_slope[[k]]$estimate*10
  period_slope[[k]]$slope_decades %>% hist()
  period_slope[[k]]$slope_decades %>% range()

}

# inspect habitat-levels estimated slopes
bind_rows(period_slope)$slope_decades %>% 
  hist
bind_rows(period_slope)$slope_decades %>% 
  range
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
means_cntrd <- means_dats %>%
  map(\(x){
  x %>%
  select(1:3) %>%
  left_join(
    x %>% filter(year==1960) %>% select(habitat, mean) %>% rename(mean.1960=mean), 'habitat'
  ) %>%
  group_by(habitat) %>%
  mutate(
    mean.change.1960 = mean-mean.1960,
    perc.change.1960 = 100*(mean.change.1960/mean.1960),
    mean.habitat = mean(mean),
    mean.change = mean-mean.habitat,
    perc.change = 100*(mean.change/mean.habitat)
  )
  })

# get the plots
p_means=list()
p_means.1960=list()
hab = c('Forest','Grassland','Scrub','Wetland')
for(k in names(means_dats)) {
  
  for (h in hab) {
  
  t4 <- period_slope[[k]] %>% filter(period=='From 1960 to 2020') %>% filter(habitat==h) %>% pull(slope_decades)
  t4 <- pals[which(vals %in% round(t4, 2))]
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  nplots = means_dats[[k]] %>% filter(habitat==h) %>% pull(n) %>% unique %>%
     prettyNum(big.mark=",",scientific=FALSE)

  ## plot means
  y_limits = dat[[k]] %>% filter(habitat==h) %>% pull(eiv) %>% quantile(c(.25,.75))
  p_means[[h]][[k]] <- means_dats[[k]] %>%
    filter(habitat==h) %>%
    ggplot(aes(year, mean)) +
    # color background
    annotate('rect', xmin=1960, xmax=2020, ymin=-Inf, ymax=Inf, alpha=.8, fill=t4) +
    # plot path
    geom_path(linewidth = 1, color='grey65') +
    # plot linear regressions per periods
    geom_smooth(data = means_dats[[k]] %>% filter(year>=1960 & year<=2020 & habitat==h), method = 'lm', color='black', se=F, linewidth=.30) +
    # "fegatelli" (cit.)
    labs(y=paste0('Mean ', k), x='Year') +
    theme_bw()+
    theme(axis.title = element_blank(),
          panel.grid = element_blank())+
    lims(y=y_limits) +
    # annotate no. plots
    annotate('text', x = 1990, y = Inf, label = paste0('n = ', nplots), vjust = 1, size= 2.75)

  ## plot relative changes over 1960 as baseline
  y_limits_centrd <- c(-.4, .4)
  p_means.1960[[h]][[k]] <- means_cntrd[[k]] %>%
    filter(habitat==h) %>%
    ggplot(aes(year, mean.change.1960)) +
    # color background
    annotate('rect', xmin=1960, xmax=2020, ymin=-Inf, ymax=Inf, alpha=.8, fill=t4) +
    # plot zero change line
    geom_segment(x = 1960, y = 0, xend = 2020, yend = 0, lty=3, color='grey30', linewidth=.3) +
    # plot path
    geom_path(linewidth = 1, color='grey65') +
    # plot linear regressions per periods
    geom_smooth(data = means_cntrd[[k]] %>% filter(year>=1960 & year<=2020 & habitat==h), method = 'lm', color='black', se=F, linewidth=.30) +
    # "fegatelli" (cit.)
    labs(y=paste0('Mean ', k), x='Year') +
    theme_bw()+
    theme(axis.title = element_blank(),
          panel.grid = element_blank())+
    # set appropriate scale 
    scale_y_continuous(
      limits = y_limits_centrd,
      labels = c( '-0.4', '-0.2',  round((means_cntrd[[k]] %>% filter(habitat==h) %>% pull(mean.1960) %>% unique()),1),  '+0.2',  '+0.4'),
      breaks = seq(-.4, .4, .2)
    ) +
    # annotate no. plots
    annotate('text', x = 1990, y = Inf, label = paste0('n = ', nplots), vjust = 1, size= 2.75)
 }
}

# prepare histogram-like column bars to put on top of trends figures
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

# get histograms
for(k in names(dat)) { 
  dk = dat_cat[[k]] 
  dk$conta <- 1:nrow(dk)
  norm.mean <- weighted.mean(dk$conta, dk$n)
  p.hists[[k]] <- dk %>%
      ggplot(aes(x = conta, y = n, fill =  mean_cat)) +
      geom_col(col=NA, position = 'dodge', width=1) +
      geom_vline(xintercept = mean(dk$conta), lty=2, col='grey50') +
      geom_vline(xintercept = norm.mean, lty=1) +
      annotate('text', x = norm.mean, y = 0, label=dat_mean.cat[[k]]) +
      scale_fill_manual(values=colpa$cols) +
      theme_bw() +
      ggtitle(k)+
      scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-3))+
      theme(legend.position = 'none', 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.title = element_blank(),
            panel.border = element_blank(), 
            axis.line.y = element_line(colour = 'black'),
            axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# combine plots into a grid
cp <- cowplot::plot_grid(
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
ggsave(paste0(pth2fig, 'average.trends.from1960.svg'), cp, width = 8, height = 7.5)

# quit
save('no')