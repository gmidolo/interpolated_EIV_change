#### Remove any plots from the main dataset from the subset where CWM was calculated #####
library(tidyverse)

username <- str_split(getwd(), '/')[[1]][3]
setwd(paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU_EIV/'))

EVA_filename = './code/data/input/EVA_01.11.2025.csv'
ReSuClean_filename  = './code/data/input/ReSurveyEU_clean_01.11.2025.csv'
EVA_ReSu_CWM_filename = './code/data/input/EVA_ReSu_CWM_01.11.2025.csv'

dat <- bind_rows(
  read_csv(EVA_filename),
  read_csv(ReSuClean_filename)
)

EVA_ReSu_CWM <- read_csv(EVA_ReSu_CWM_filename) 
nrow(EVA_ReSu_CWM)

EVA_ReSu_CWM <- EVA_ReSu_CWM %>%
  semi_join(dat, 'plot_id')
nrow(EVA_ReSu_CWM)

dat_CWM <- dat %>% 
  semi_join(EVA_ReSu_CWM, 'plot_id') %>%
  left_join(
    EVA_ReSu_CWM %>% select(plot_id, contains('cwm.EIV'))
  )
names(dat_CWM)
dat_CWM <- dat_CWM %>% select(1:18, contains('cwm.'), everything()) # reorder columns

glimpse(dat_CWM)

# just check resurveys
dat_CWM_RSRV <- dat_CWM %>% filter(!is.na(resurv_id))
dat_CWM_RSRV %>% group_by(resurv_id) %>% summarise(n=n()) %>% pull(n) %>% range()

## Export
write_csv(dat_CWM, file = EVA_ReSu_CWM_filename)



## test of correlation between weighted vs. unweighted methods
nrow(dat_CWM)
# correlations across a total of 572,241 plots (as for 29.11.25)
eiv_codes <- c('M', 'N', 'R', 'L', 'T')
plot_titles <- c('Moisture', 'Nitrogen', 'Reaction', 'Light', 'Temperature')
plot_eiv_corr <- function(code, title_text, data) {
  x_var <- paste0('cwm.EIV_', code)
  y_var <- paste0('cm.EIV_', code)
  di <- data[,c(x_var,y_var)] %>% 
    drop_na()
  p <- ggplot(di, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_hex(bins=40) +
    coord_fixed() +
    geom_abline(intercept = 0, slope = 1, lwd=.75, col='grey40', lty=3, alpha=1) +
    geom_smooth(col='brown3', alpha=.5, lwd = 0.75) +
    scale_fill_continuous(low = 'gray90', high = 'midnightblue', trans = scales::log10_trans()) +
    theme_bw() +
    ggtitle(title_text) +
    labs(fill='no. plots', x='EIV Community Weighted Mean', y='EIV Community Mean') +
    ggpubr::stat_cor(method = 'pearson', label.x = 1, label.y = 9) +
    annotate("text", x = 3.5, y = 8.5, label = paste0('no. plots = ', prettyNum(nrow(di),big.mark=','))) +
    scale_x_continuous(breaks = c(0,2.5,5,7.5,10), limits = c(0, 10)) +
    scale_y_continuous(breaks = c(0,2.5,5,7.5,10), limits = c(0, 10)) +
    theme(legend.position = 'bottom')
  
  return(p)
}
eiv_corr_plots <- mapply(
  FUN = plot_eiv_corr,
  code = eiv_codes,
  title_text = plot_titles,
  MoreArgs = list(data = EVA_ReSu_CWM), # Pass CWM_eive as a fixed argument
  SIMPLIFY = FALSE
)
names(eiv_corr_plots) <- plot_titles

for (i in names(eiv_corr_plots)) {
  ggsave(paste0('./fig/CWM_CM_corr_',i,'.pdf'), eiv_corr_plots[[i]], width = 4.2, height = 4.6)
}