library(tidyverse)

username <- Sys.info()['user']
setwd(paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU_EIV/'))

## load needed data in EVA, ReSurveyEurope and EIVE species
EVA = './code/data/input/EVA_01.11.2025.csv' %>% read_csv()
ReSuClean  = './code/data/input/ReSurveyEU_clean_01.11.2025.csv' %>% read_csv()
EIVE_EVA_merge <- './code/data/eiv.data.nomenclature/EIVE_EVA_merge_01.11.2025.csv' %>% read_csv()
drevojan <- paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code/data/eiv.data.nomenclature/life_form.xlsx') %>%
  readxl::read_excel() %>%
  rename(matched_concept = FloraVeg.Taxon)

## Load core EVA species data and process it
species_data <- read_delim(
  paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU/data/eva/222_PlantDiversity20241017_notJUICE/222_PlantDiversity20241017_notJUICE_species.csv'),
  col_select = c(1,3,6,9,10), show_col_types = F
) %>%
  bind_rows( ## add GVRD species data
    read_delim(
      paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU/data/eva/222_PlantDIversity20241025_GVRD/222_GVRD_20241025_notJUICE_species.csv'),
      col_select = c(1,3,6,9,10), show_col_types = F
    ) 
  ) %>%
  bind_rows( ## add UKFloodplainMeadows species data
    read_delim(
      paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU/data/eva/222_UKFloodplainMeadows20241031/222_UKFloodplainMeadows20241031_notJUICE_species.csv'),
      col_select = c(1,3,6,9,10), show_col_types = F
    ) 
  ) %>%
  setNames(c('plot_id','taxongroup', 'matched_concept', 'layer', 'cover')) %>% # rename columns
  semi_join(bind_rows(EVA,ReSuClean), 'plot_id') # filter out only the plots we selected this far

# add artificial cover (1) to Denmark Nature data
DKnature_plots <- EVA %>% filter(dataset=='Denmark Naturdata') %>% select(plot_id)
species_data$cover <- ifelse(species_data$plot_id %in% DKnature_plots$plot_id, 1, species_data$cover)

#only include species with cover >0
species_data <- species_data %>% filter(cover>0)

## check taxa unknown in EVA: are they vascular plants?
unk_taxa_freq <- species_data %>% 
  filter(taxongroup == 'Unknown') %>% 
  select(plot_id, matched_concept) %>% ungroup() %>% unique() %>%
  group_by(matched_concept) %>%
  summarise(n = n()) # check frequency of Unknown taxa in the dataset
unk_taxa <- species_data %>% 
  filter(taxongroup == 'Unknown') %>% 
  select(matched_concept) %>% 
  unique() 
unk_genera <- unk_taxa %>% 
  mutate(matched_concept= str_remove(str_remove(matched_concept,'Cf. '),'cf. ') %>% str_remove('cf ')) %>% 
  separate(matched_concept, into=c('genus'), sep=' ') %>%
  mutate(genus=str_remove_all(genus, '\\*'))
#remotes::install_github('helixcn/plantlist') # install plantlist R package
check_unk_genera <- plantlist::TPL(unk_genera$genus) # search the genus for higher plants under modern classification systems
check_unk_genera <- unique(check_unk_genera) %>% select(YOUR_SEARCH,GROUP) %>% setNames(c('genus', 'GROUP'))
check_unk_genera$GROUP <- ifelse(is.na(check_unk_genera$GROUP), 'NotFound', check_unk_genera$GROUP)
unk_genera <- unk_genera %>% left_join(check_unk_genera, 'genus')
unk_genera$keep_unk_genera  <- unk_genera$GROUP %in% (c('Angiosperms', 'Ferns and lycophytes', 'Gymnosperms'))
unk_taxa_2keep <- unk_taxa[unk_genera$keep_unk_genera, ]

# list of vascular taxa (will be used to filter the vascular plant upon which to calculate species richness)
vascular_taxa <- species_data %>%
  filter(taxongroup == 'Vascular plant') %>% 
  select(matched_concept) %>%
  unique() %>%
  bind_rows(unk_taxa_2keep) %>% # bind checked names
  arrange(matched_concept)

# finally, select species that are vascular plants only
species_data <- species_data %>%
  semi_join(vascular_taxa, 'matched_concept')

## 1. No trees ####
drevojan_treeshrub <- drevojan %>% filter(Tree==1 | Shrub==1) %>% select(matched_concept)
EIVE_EVA_merge_treeshrub <- EIVE_EVA_merge %>% semi_join(drevojan_treeshrub)
nrow(EIVE_EVA_merge_treeshrub) # no. trees excluded: 845
EIVE_EVA_merge_NO_treeshrub <- EIVE_EVA_merge %>% anti_join(EIVE_EVA_merge_treeshrub, 'matched_concept')

species_data_eive <- species_data %>%
  select(plot_id, matched_concept) %>%
  ungroup() %>%
  unique() %>%
  anti_join(EIVE_EVA_merge_treeshrub) %>%
  left_join(EIVE_EVA_merge_NO_treeshrub)

glimpse(species_data_eive)

# calc cm (community means)
n_na <- \(x){sum(rep(1, length(x))[!is.na(x)])}
mean_na <- \(x){mean(x, na.rm = T)}
sd_na <- \(x){sd(x, na.rm = T)}

st=Sys.time() # approx. 1.5 mins
CM_eive_NOTREES <- species_data_eive %>%
  group_by(plot_id) %>%
  summarise(
    n=n(),
    n.EIV_M=n_na(EIV_M)/n,
    n.EIV_N=n_na(EIV_N)/n,
    n.EIV_R=n_na(EIV_R)/n,
    n.EIV_L=n_na(EIV_L)/n,
    n.EIV_T=n_na(EIV_T)/n,
    cm.EIV_M=mean_na(EIV_M),
    cm.EIV_N=mean_na(EIV_N),
    cm.EIV_R=mean_na(EIV_R),
    cm.EIV_L=mean_na(EIV_L),
    cm.EIV_T=mean_na(EIV_T),
    sd.EIV_M=sd_na(EIV_M),
    sd.EIV_N=sd_na(EIV_N),
    sd.EIV_R=sd_na(EIV_R),
    sd.EIV_L=sd_na(EIV_L),
    sd.EIV_T=sd_na(EIV_T)
  ) %>%
  ungroup()
Sys.time()-st

par(mfrow=c(2,3))
hist(CM_eive_NOTREES$cm.EIV_M)
hist(CM_eive_NOTREES$cm.EIV_N)
hist(CM_eive_NOTREES$cm.EIV_R)
hist(CM_eive_NOTREES$cm.EIV_L)
hist(CM_eive_NOTREES$cm.EIV_T)

table(CM_eive_NOTREES$n.EIV_M >= 0.8) %>% prop.table()
table(CM_eive_NOTREES$n.EIV_N >= 0.8) %>% prop.table()
table(CM_eive_NOTREES$n.EIV_R >= 0.8) %>% prop.table()
table(CM_eive_NOTREES$n.EIV_L >= 0.8) %>% prop.table()
table(CM_eive_NOTREES$n.EIV_T >= 0.8) %>% prop.table()


EVA_ReSuEU <- bind_rows(EVA, ReSuClean) %>% arrange(plot_id)

EVA_ReSuEU_NOTREES <- CM_eive_NOTREES %>%
  left_join(
    EVA_ReSuEU %>%
      semi_join(CM_eive_NOTREES, 'plot_id') %>%
      select(-c(n:sd.EIV_T))
  ) %>%
  select(
    database, plot_id, dataset, ESy, habitat, lon, lat, x, y, elev, year, plot_size, everything()
  ) %>% arrange(plot_id)

EVA_ReSuEU_semi = EVA_ReSuEU %>% semi_join(EVA_ReSuEU_NOTREES, 'plot_id')
table(EVA_ReSuEU_semi$plot_id == EVA_ReSuEU_NOTREES$plot_id)

ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

di = list()
for (i in ind.names$eiv_name_raw) {
  di[[i]] <- cbind(EVA_ReSuEU_NOTREES[[paste0('cm.',i)]], EVA_ReSuEU_semi[[paste0('cm.',i)]]) %>%
    as.data.frame() 
  di[[i]]$h = EVA_ReSuEU_semi$habitat
  di[[i]]$plot_id = EVA_ReSuEU_semi$plot_id
  di[[i]] = di[[i]] %>% drop_na()
}
di <- di %>% bind_rows(.id = 'eiv_name_raw') %>% left_join(ind.names)
di$eiv_name <- factor(di$eiv_name, c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))
di_counts <- di %>%
  group_by(h, eiv_name) %>%
  summarise(n_obs = n()) %>%
  ungroup()

pi <- di %>%
  ggplot(aes(V1,V2)) +
  geom_hex(bins=40) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1, lwd=.75, col='grey40', lty=3, alpha=1) +
  geom_smooth(col='brown3', alpha=.5, lwd = 0.75) +
  scale_fill_continuous(low = 'gray90', high = 'midnightblue', trans = scales::log10_trans()) +
  theme_bw() +
  labs(fill='no. plots',
       y= bquote(CM[EIV]~" - all species included"),
       x = bquote(CM[EIV]~" - tree and shrub species excluded")
  ) +
  ggpubr::stat_cor(method = 'pearson', label.x = 0.1, label.y = 9, size = 2.75) +
  geom_text(
    data = di_counts, 
    aes(label = paste0('no. plots = ', prettyNum(n_obs, big.mark=','))), x = 3.5, y = 8.5,        
    inherit.aes = FALSE, size = 2.75
  ) +
  scale_x_continuous(breaks = c(0,2.5,5,7.5,10), limits = c(0, 10)) +
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10), limits = c(0, 10)) +
  theme(legend.position = 'bottom') +
  facet_grid(h~eiv_name)

ggsave(paste0('./fig/notrees_CM_corr.pdf'), pi, width = 10, height = 8)

write_csv(EVA_ReSuEU_NOTREES, file = './code/data/input/EVA_ReSu_NOTREES_01.11.2025.csv')


# ## 2. Non-rare species only ####
rarity.threshold <- c(1500)

EIVE_EVA_NONrare <- EIVE_EVA_merge %>% filter(freq >= rarity.threshold)

species_data_eive <- species_data %>%
  select(plot_id, matched_concept) %>%
  ungroup() %>%
  unique() %>%
  left_join(EIVE_EVA_NONrare)
glimpse(species_data_eive)

# calc cm (community means)
n_na <- \(x){sum(rep(1, length(x))[!is.na(x)])}
mean_na <- \(x){mean(x, na.rm = T)}
sd_na <- \(x){sd(x, na.rm = T)}

st=Sys.time() # approx. 1.5 mins
CM_eive_NONrare <- species_data_eive %>%
  group_by(plot_id) %>%
  summarise(
    n=n(),
    n.EIV_M=n_na(EIV_M)/n,
    n.EIV_N=n_na(EIV_N)/n,
    n.EIV_R=n_na(EIV_R)/n,
    n.EIV_L=n_na(EIV_L)/n,
    n.EIV_T=n_na(EIV_T)/n,
    cm.EIV_M=mean_na(EIV_M),
    cm.EIV_N=mean_na(EIV_N),
    cm.EIV_R=mean_na(EIV_R),
    cm.EIV_L=mean_na(EIV_L),
    cm.EIV_T=mean_na(EIV_T),
    sd.EIV_M=sd_na(EIV_M),
    sd.EIV_N=sd_na(EIV_N),
    sd.EIV_R=sd_na(EIV_R),
    sd.EIV_L=sd_na(EIV_L),
    sd.EIV_T=sd_na(EIV_T)
  ) %>%
  ungroup()
Sys.time()-st

glimpse(CM_eive_NONrare)
CM_eive_NONrare_cmonly <- CM_eive_NONrare %>% select(plot_id, contains('cm.')) %>%
  gather('k','vrare',cm.EIV_M:cm.EIV_T) 

EVA_orig <- bind_rows(EVA,ReSuClean) %>% select(plot_id, contains('cm.')) %>%
  gather('k','vorig',cm.EIV_M:cm.EIV_T) 

ind.names <- data.frame(
  eiv_name_raw = c('EIV_L','EIV_T','EIV_M','EIV_N','EIV_R'),
  eiv_name = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

d2p <- EVA_orig %>% left_join(CM_eive_NONrare_cmonly) %>%
  mutate(eiv_name_raw = str_remove(k, 'cm.')) %>%
  left_join(ind.names)

d2p <- d2p %>% drop_na()

d2p$eiv_name <- factor(d2p$eiv_name, c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction'))

d2p_count <- d2p %>% group_by(eiv_name) %>% summarise(n = prettyNum(n(),big.mark=','))
d2p_count

p = ggplot(d2p, aes(x=vrare,y=vorig)) +
  geom_hex(bins=40) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1, lwd=.75, col='grey40', lty=3, alpha=1) +
  # geom_smooth(col='brown3', alpha=.5, lwd = 0.75) +
  scale_fill_continuous(low = 'gray90', high = 'midnightblue', trans = scales::log10_trans()) +
  theme_bw() +
  labs(fill='no. plots') +
  ggpubr::stat_cor(method = 'pearson', label.x = 1, label.y = 9, size = 2.75) +
  labs(y=expression(CM[EIV]~'- without "rare" species (1,258 species)'), x=expression(CM[EIV]~'- all species (13,874 species)')) +
  #annotate("text", x = 3.5, y = 8.5, label = paste0('no. plots = ', prettyNum(nrow(di),big.mark=','))) +
  scale_x_continuous(breaks = c(0,2.5,5,7.5,10), limits = c(0, 10)) +
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10), limits = c(0, 10)) +
  facet_wrap(~eiv_name)
p
  
ggsave("C:/Users/gabri/OneDrive - CZU v Praze/czu/intrplEU_EIV/ms_SciAdv/resub/Rplot_rare_species_missing.png",p, width = 8, height = 4, dpi = 500)
