#### 1. Load data for duplication check #####
library(tidyverse);library(sf)

username <- str_split(getwd(), '/')[[1]][3]
setwd(paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU_EIV/'))

eive_eva_merge_filename = './code/data/eiv.data.nomenclature/EIVE_EVA_merge_01.11.2025.csv'
EVA_filename = './code/data/input/EVA_01.11.2025.csv'
ReSuEU_filename  = './code/data/input/ReSurveyEU_01.11.2025.csv'

## Load cleaned EVA data
EVA <- EVA_filename %>%
  read_csv(show_col_types = F)

## Load cleaned ReSuEU
ReSuEU <- ReSuEU_filename %>%
  read_csv(show_col_types = F)

## Load core EVA species data
species_data <- read_delim(
  paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU/','data/eva/222_PlantDiversity20241017_notJUICE/222_PlantDiversity20241017_notJUICE_species.csv'),
  col_select = c(1,3,6,9,10), show_col_types = F
) %>%
  bind_rows( ## add GVRD species data
    read_delim(
      paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU/','data/eva/222_PlantDiversity20241017_notJUICE/222_PlantDiversity20241017_notJUICE_species.csv'),
      col_select = c(1,3,6,9,10), show_col_types = F
    ) 
  ) %>%
  bind_rows( ## add UKFloodplainMeadows species data
    read_delim(
      paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU/','data/eva/222_PlantDiversity20241017_notJUICE/222_PlantDiversity20241017_notJUICE_species.csv'),
      col_select = c(1,3,6,9,10), show_col_types = F
    ) 
  ) %>%
  setNames(c('plot_id','taxongroup', 'matched_concept', 'layer', 'cover'))  %>% # rename columns
  semi_join(data.frame(plot_id = c(EVA$plot_id, ReSuEU$plot_id)), 'plot_id') # filter out only the plots we selected this far

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

## list of vascular taxa (will be used to filter the vascular plant upon which to calculate species richness)
vascular_taxa <- species_data %>%
  filter(taxongroup == 'Vascular plant') %>% 
  select(matched_concept) %>%
  unique() %>%
  bind_rows(unk_taxa_2keep) %>% # bind checked names
  arrange(matched_concept)

## Finalize species data (keep vascular taxa only)
species_data <-  species_data %>%
  semi_join(vascular_taxa, 'matched_concept')

species_data$cover <- ifelse(species_data$cover == 0, 100, species_data$cover) # put 100 cover in plots that are missing cover values

# Standardize cover via sum
species_data <- species_data %>%
  group_by(plot_id, matched_concept) %>%
  summarise(sum.cover = sum(cover)) %>%
  ungroup()

species_list <- species_data %>% 
  filter(!str_detect(matched_concept, ' species')) %>%
  group_by(matched_concept) %>% 
  summarise(frq=n()) %>% 
  ungroup()

# 5. check selected EIVE species after merging from Irena ####
##!!!!Achtung!!!!##
#PLAY IT AT THE END OF DATA CLEANING (after running duplicate search and ReSurvey final clean) FOR ACCURATE ESTIMATES#
eive_eva_merge <- eive_eva_merge_filename %>%
  read_csv(show_col_types = F)

eive_eva_merge %>% semi_join(species_list) %>% nrow() # how many species did we actually use? 13700 # values updated on 01 Nov 2025
species_list %>% anti_join(eive_eva_merge) %>% nrow() # how many species we coundnt merge? 2544
species_list %>% nrow() # full vascular species available in eva (sum of the previous two): 16244
species_list %>% anti_join(eive_eva_merge) %>% pull(frq) %>% range %>% round()
species_list %>% anti_join(eive_eva_merge) %>% pull(frq) %>% mean %>% round()
species_list %>% anti_join(eive_eva_merge) %>% pull(frq) %>% sd %>% round()

#2. Duplication detection WITHIN: EVA #####

## Value in meters to aggregate coordinates
fact = 1000 # 1 km

## Aggregate in space
EVA_aggr <- EVA %>%
  mutate(x_aggr = round(x/fact, 0)*fact,
         y_aggr = round(y/fact, 0)*fact) %>%
  group_by(x_aggr,y_aggr) %>%
  mutate(ID_aggr=cur_group_id(), .before = plot_id) %>%
  ungroup()

length(unique(EVA_aggr$ID_aggr))

## Flag records within each spatial aggregate
EVA_flag <- EVA_aggr %>%
  group_by(n, year, plot_size, ESy, ID_aggr) %>% # define here the grouping parameters
  mutate(flag = n() > 1) %>%
  mutate(ID_flag=cur_group_id(), .before=plot_id) %>%
  ungroup()

EVA_flag$flag %>% table()

EVA_flag <- EVA_flag %>%
  filter(flag) 

## Check for absolute replication (exact same species composition and cover)
st = Sys.time()
EVA_flag_compo <- EVA_flag %>%
  select(plot_id, ID_flag) %>%
  left_join(species_data, 'plot_id') %>%
  spread(matched_concept, sum.cover, fill=0) %>%
  group_by(across(c(-plot_id))) %>%
  mutate(totrep = n() > 1) %>%
  mutate(ID_totrep = cur_group_id(), .before=plot_id)  %>%
  ungroup() %>%
  select(plot_id, totrep, ID_totrep)
Sys.time()-st # ~ 1.2 mins!!!

# EVA_duplicated_for_Ilona <- 
#   read_delim('./data/eva/222_PlantDiversity20241017_notJUICE/222_PlantDiversity20241017_notJUICE_header.csv') %>%
#   mutate(Altitude=as.numeric(Altitude), `Slope (°)`=as.numeric(`Slope (°)`)) %>%
#   bind_rows(read_delim('./data/eva/222_PlantDIversity20241025_GVRD/222_GVRD_20241025_notJUICE_header.csv')) %>%
#   semi_join(EVA_flag_compo %>% filter(totrep), by = join_by('PlotObservationID' == 'plot_id')) %>%
#   left_join(EVA_flag_compo %>% filter(totrep) %>% select(1,3) %>% rename(Duplication_ID=ID_totrep), by = join_by('PlotObservationID' == 'plot_id')) %>%
#   select(Duplication_ID, everything()) %>%
#   arrange(Duplication_ID)
# 
# write_delim(EVA_duplicated_for_Ilona, file = './data/eva/EVA_flagged_duplicates.txt', delim = '\t', )

## Randomply pick (presumed) duplicated records to be excluded
set.seed(1975)
EVA_flag_keep <- EVA_flag_compo %>% 
  group_by(ID_totrep) %>% 
  sample_n(1) %>% 
  ungroup()
EVA_flag_remove <- EVA_flag_compo %>% 
  anti_join(EVA_flag_keep)


## Map plots to be removed
EU <- read_sf(list.files(paste0('C:/Users/',username,'/OneDrive - CZU v Praze/brno/brno_postdoc/maps/euro+med_map/Final.Map'), pattern = 'shp',full.names = T)) %>%
  semi_join(data.frame(name =
                         c('Albania', 'Austria', 'Baleares', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria',
                           'Corsica', 'Crete', 'Croatia', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany',
                           'Greece', 'Hungary', 'Ireland', 'Italy', 'Kaliningrad Region', 'Kosovo', 'Latvia', 'Liechtenstein',
                           'Lithuania', 'Luxembourg', 'Malta', 'Moldova', 'Montenegro', 'Netherlands', 'North Macedonia',
                           'Norway', 'Poland', 'Portugal', 'Romania', 'Sardinia', 'Serbia', 'Sicily',
                           'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine', 'United Kingdom'
                         )), 'name') %>%
  st_transform(crs = 25832) %>%
  st_simplify(dTolerance = 1000)
p <- ggplot() + 
  geom_sf(data = EU, fill='white') +
  geom_hex(data=EVA_aggr %>%
             semi_join(EVA_flag_remove), aes(x_aggr, y_aggr), bins = 90, alpha = 0.8) +
  scale_fill_viridis_c(trans='log10', option='C') +
  theme_bw() +
  theme(axis.title = element_blank()) +
  labs(fill = 'N. of duplicated \nplots in EVA to be removed')
p

# Finalize EVA!
EVA <- EVA %>%
  anti_join(EVA_flag_remove, 'plot_id')


#3. Duplication detection ACROSS: EVA and ReSurveyEU ####
#value in meters to aggregate coordinates
fact = 1000  

## Bind databases rows
EVAReSuEU <- bind_rows(EVA, ReSuEU)

## Aggregate in space
EVAReSuEU_aggr <- EVAReSuEU %>%
  mutate(x_aggr = round(x/fact, 0)*fact,
         y_aggr = round(y/fact, 0)*fact) %>%
  group_by(x_aggr,y_aggr) %>%
  mutate(ID_aggr=cur_group_id(), .before = plot_id) %>%
  ungroup()

length(unique(EVAReSuEU_aggr$ID_aggr)) # number of "cells"

## Flag records within each spatial aggregate
EVAReSuEU_flag <- EVAReSuEU_aggr %>%
  group_by(n, year, plot_size, ESy, ID_aggr) %>% # define here the grouping parameters
  mutate(flag = n() > 1) %>%
  mutate(ID_flag=cur_group_id(), .before=plot_id) %>%
  ungroup()

EVAReSuEU_flag$flag %>% table()

EVAReSuEU_flag <- EVAReSuEU_flag %>%
  filter(flag) 

## Filter flagged plots with same ID, but which are detected across the two databases
EVAReSuEU_flag <- EVAReSuEU_flag %>%
  ungroup() %>%
  group_by(ID_flag) %>%
  mutate(ID_bothDBs=cur_group_id(), .before=plot_id) %>%
  mutate(countDBs = n_distinct(database), .before=plot_id) %>%
  mutate(bothDBs = countDBs > 1, .before=plot_id) %>%
  ungroup()

EVAReSuEU_flag$bothDBs %>% table()
EVAReSuEU_flag$countDBs %>% table()

EVAReSuEU_flag <- EVAReSuEU_flag %>%
  filter(bothDBs) %>%
  arrange(ID_bothDBs)

# EVAReSuEU_flag %>% group_by(ID_bothDBs) %>% summarise(n=n()) %>% View()


## Check for absolute replication (exact same species composition and cover)
st = Sys.time()
EVAReSuEU_flag_compo <- EVAReSuEU_flag %>%
  select(plot_id, ID_bothDBs) %>%
  left_join(species_data, 'plot_id') %>%
  spread(matched_concept, sum.cover, fill=0) %>%
  group_by(across(c(-plot_id))) %>%
  mutate(totrep = n() > 1, .before=plot_id) %>%
  mutate(ID_totrep = cur_group_id(), .before=plot_id)  %>%
  ungroup() %>%
  select(plot_id, totrep, ID_totrep)
Sys.time()-st

## Identify duplicated plots in ReSurveyEU that are identical to plots already present in EVA. They will be discarded from ReSurveyEU
to_exclude_in_ReSurveyEU = EVAReSuEU_flag %>%
  semi_join(EVAReSuEU_flag_compo %>% filter(totrep)) %>%
  filter(database == 'ReSurveyEU')

ReSuEU <- ReSuEU %>%
  anti_join(to_exclude_in_ReSurveyEU, 'plot_id')

#4. Overwrite / export the data ####

write_csv(EVA, file = EVA_filename)
write_csv(ReSuEU, file = ReSuEU_filename)