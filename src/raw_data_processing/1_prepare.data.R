# 1. load pkgs####
library(tidyverse);library(sf);library(terra)

username <- str_split(getwd(), '/')[[1]][3]
setwd(paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU/'))

date.update='01.11.2025'

# 2. Load header data####

renames <- c('plot_id', 'dataset', 'lon', 'lat', 'country', 'plot_size', 'ESy', 'full_recording_date', 'location_uncertainty', 'for_EVA_y_n', 'ReSurvey_y_n', 'ReSurvey_type', 'ReSur_site', 'ReSur_plot', 'ReSur_obs', 'manipulate_y_n', 'ReSur_time','ReSur_dupl','manipulation_type')

## read EVA + ReSurvey 222 core data
EVA_ReSu <- 
  './data/eva/222_PlantDiversity20241017_notJUICE/222_PlantDiversity20241017_notJUICE_header.csv' %>%
  read_delim(col_select = c('PlotObservationID','Dataset', 'Longitude', 'Latitude', 'Country', 'Relevé area (m²)', 'Expert System','Date of recording', 'Location uncertainty (m)', 'For EVA (Y/N)', 'ReSurvey plot (Y/N)',
                            'RS_PROJTYP', 'ReSurvey site', 'ReSurvey plot', 'ReSurvey observation', 'Manipulate (y/n)','RS_TIME','RS_DUPL','Type of manipulation')
             ) 
## add GVRD-Germany data
EVA_ReSu <- bind_rows(EVA_ReSu,
                    './data/eva/222_PlantDIversity20241025_GVRD/222_GVRD_20241025_notJUICE_header.csv' %>%
                      read_delim(col_select = colnames(EVA_ReSu), show_col_types = F)
                      )

## add UKFloodplainMeadows data
EVA_ReSu <- bind_rows(EVA_ReSu,
                      './data/eva/222_UKFloodplainMeadows20241031/222_UKFloodplainMeadows20241031_notJUICE_header.csv' %>%
                        read_delim(col_select = colnames(EVA_ReSu), show_col_types = F) %>%
                        mutate(`ReSurvey site` = as.character(`ReSurvey site`))
                      )

## rename colnames
EVA_ReSu <- setNames(EVA_ReSu, renames)

## subset Denmark Naturdata metadata 
EVA_ReSu_DN <- EVA_ReSu %>%
  filter(str_detect(dataset,'Denmark Natur'))

## clean EVA_ReSu from old Denmark data and bind the new version to it
# first, read Denmark Naturdata data with habitat...
rawDenNat <- './data/eva/Naturdata_Eunis-RSy.csv' %>%
  read_delim(delim=';', show_col_types = F,
             col_select = c('PlotObservationID','HABITAT','Eunis-ESy', 'EUNIS-ESy level 1')
             ) 
rawDenNat <- setNames(rawDenNat, c('plot_id','DenNatdat_habitat','ESy_max', 'ESy1'))
rawDenNat$ESy <- ifelse(is.na(rawDenNat$ESy_max), rawDenNat$ESy1, rawDenNat$ESy_max)

# ...then,isolate non-danish data...
EVA_ReSu <- EVA_ReSu %>%
  anti_join(EVA_ReSu_DN, 'plot_id')

# ...then, correct habitat errors in non-danish data (Ilona correction from Oct 2025)...
correct_header <- './data/eva/222_OnlyEUNISclass20251001/222_OnlyEUNISclass20251001.txt' %>%
  read_delim(show_col_types = F) %>%
  setNames(c('plot_id', 'ESy_Ilona.October.2025')) %>%
  anti_join(EVA_ReSu_DN, 'plot_id')
table(EVA_ReSu$plot_id %in% correct_header$plot_id)
EVA_ReSu <- EVA_ReSu %>%
  left_join(
    correct_header, 'plot_id'
  )
prop.table(table(EVA_ReSu$ESy_Ilona.October.2025 == EVA_ReSu$ESy))
EVA_ReSu <- EVA_ReSu %>%
  mutate(ESy = ESy_Ilona.October.2025) %>%
  select(-ESy_Ilona.October.2025)

# ...finally, bind the non-danish data to danish ones
EVA_ReSu_DN <- EVA_ReSu_DN %>% 
  select(-ESy) %>% 
  left_join(rawDenNat %>% select(plot_id, ESy), 'plot_id')
EVA_ReSu <- EVA_ReSu %>%
  bind_rows(EVA_ReSu_DN)

## arrange by plot id
EVA_ReSu <- EVA_ReSu %>%
  arrange(plot_id)


# 3. Apply standard EVA filters (habitat, plot size, and year) #### 

## load habitat conversion table
habitat_conversion <- 
  paste0('./data/eva/habitat_conversion_table.csv') %>%
  read_csv(show_col_types = F) %>%
  arrange(habitat, ESy) %>%
  semi_join(data.frame(habitat=c('forest','grassland','scrub','wetland'))) %>%
  mutate(habitat = str_to_title(habitat)) %>%
  select(ESy, habitat)
head(habitat_conversion)

## filter habitat needed and categorize them into 'forest','grassland','scrub',or 'wetland'
EVA_ReSu <- EVA_ReSu %>% 
  semi_join(habitat_conversion, 'ESy') %>%
  left_join(habitat_conversion, 'ESy')

EVA_ReSu %>% # explore distribution of the data (plot sizes)
  ggplot(aes(x=plot_size)) +
  geom_histogram() +
  scale_x_continuous(trans = 'log10', breaks = c(1, 10 ,100, 1000))+
  facet_wrap(~habitat, scales = 'free') +
  theme_bw()

## Define filtering table
plot_size_table <- data.frame(
  habitat = c('Forest', 'Grassland', 'Scrub', 'Wetland'),
  min_size = c(100, 1, 1, 1),
  max_size = c(1000, 100, 100, 100)
)

## Apply dedicated filter
keep_size <- list()
for (i in plot_size_table$habitat) {
  keep_size[[i]] <- EVA_ReSu %>%
    filter(habitat == i) %>%
    mutate(keep = (plot_size >= plot_size_table[which(plot_size_table$habitat %in% i), 'min_size'] &
                  plot_size <= plot_size_table[which(plot_size_table$habitat %in% i), 'max_size'])
           )

  message(i)
  print(round(prop.table(table(keep_size[[i]]$keep))*100, 2))
  
  keep_size[[i]] <- keep_size[[i]] %>% filter(keep) %>% select(plot_id) 
}
EVA_ReSu <- EVA_ReSu %>%
  semi_join(bind_rows(keep_size))

## Visually check again the distribution of plots
EVA_ReSu %>% # re-explore distribution of the data (plot sizes)
  ggplot(aes(x=plot_size)) +
  geom_histogram(bins = 10) +
  scale_x_continuous(trans = 'log10', breaks = c(1, 10 ,100, 1000))+
  facet_wrap(~habitat, scales = 'free') +
  theme_bw() + labs(x='Plot size (squared m)', y='Plot count')

## Extract year
EVA_ReSu$year <- lubridate::year(as.Date(EVA_ReSu$full_recording_date, format = "%d.%m.%Y"))
prop.table(table(!is.na(EVA_ReSu$year)))# proportion of plots that have a date
# barplot(table(lubridate::wday(as.Date(EVA_ReSu$full_recording_date, format = "%d.%m.%Y"), label=T))) # Easter egg: which day of the week is the most sampled?
# barplot(table(lubridate::month(as.Date(EVA_ReSu$full_recording_date, format = "%d.%m.%Y"), label=T))) # Most sampled months, here you can see that we have several outliers in January

# Retain only plots with sampling year
EVA_ReSu <- EVA_ReSu %>% 
  filter(!is.na(year))

# Retain only plots with sampling year >= 1945
EVA_ReSu <- EVA_ReSu %>% 
  filter(year >= 1945)

EVA_ReSu %>%
  ggplot(aes(x=year)) +
  geom_histogram() + theme_bw()

## Inspect location uncertainty (we keep everything below 1 km)
EVA_ReSu$location_uncertainty %>% hist(main = 'Distribution of location uncertainty', xlab = 'Uncertainty (m)')
table(EVA_ReSu$location_uncertainty <= 1000 | is.na(EVA_ReSu$location_uncertainty)) # plots with location unc. less than 1 km; or NA
EVA_ReSu <- EVA_ReSu %>%
  filter(
    location_uncertainty <= 1000 | is.na(location_uncertainty)
  )

## all plots must have lon and lat
EVA_ReSu <- EVA_ReSu %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(lat))


# 4. Clip out geographical outliers ####

## Define projection
prj = 25832 

## Define EU study area
EU <- read_sf(list.files(paste0('C:/Users/',username,'/OneDrive - CZU v Praze/brno/brno_postdoc/maps/euro+med_map/Final.Map'), pattern = 'shp',full.names = T)) %>%
  semi_join(data.frame(name =
                         c('Albania', 'Austria', 'Baleares', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria',
                           'Corsica', 'Crete', 'Croatia', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany',
                           'Greece', 'Hungary', 'Ireland', 'Italy', 'Kaliningrad Region', 'Kosovo', 'Latvia', 'Liechtenstein',
                           'Lithuania', 'Luxembourg', 'Malta', 'Moldova', 'Montenegro', 'Netherlands', 'North Macedonia',
                           'Norway', 'Poland', 'Portugal', 'Romania', 'Sardinia', 'Serbia', 'Sicily',
                           'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine', 'United Kingdom'
                         )), 'name') %>%
  st_transform(crs = prj) %>%
  st_simplify(dTolerance = 1000)

## Define various boundaries
EU_51 <-  EU %>% st_buffer(51*1000) # 51 km buffer
EU_51 <- st_as_sf(st_union(EU_51$geometry)) # apply union
EU_1 <-  EU %>% st_buffer(1*1000) # 1 km buffer
EU_1 <- st_as_sf(st_union(EU_1$geometry)) # apply union

## Spatialize EVA_ReSu
EVA_ReSu_spat <- EVA_ReSu %>%
  st_as_sf(coords=c('lon','lat'), crs='WGS84', remove = F) %>%
  st_transform(crs = prj)

st = Sys.time()
out_of_bounds_51km <- st_disjoint(EVA_ReSu_spat, EU_51, sparse = F)
out_of_bounds_51km <- EVA_ReSu_spat[out_of_bounds_51km[,1], ]
Sys.time()-st

nrow(out_of_bounds_51km)
plot(out_of_bounds_51km %>% filter(lat>0) %>% pull(geometry), col='red', pch=16)
plot(EU, add=T, col=NA)

EVA_ReSu_spat <- EVA_ReSu_spat %>% 
  anti_join(out_of_bounds_51km %>% st_drop_geometry(), 'plot_id')

st = Sys.time()
out_of_bounds_1km <- st_disjoint(EVA_ReSu_spat, EU_1, sparse = F)
out_of_bounds_1km <- EVA_ReSu_spat[out_of_bounds_1km[,1], ]
Sys.time()-st

out_Ukraine <- out_of_bounds_1km %>%
  filter(country == 'Ukraine') %>%
  filter(lat <= 44.6)

out_Turkey <- out_of_bounds_1km %>%
  filter(country == 'Turkey')

out_Bulgaria <- out_of_bounds_1km %>%
  filter(country == 'Bulgaria')

additional2remove <- data.frame(plot_id = c(999121,999139,901878,1124870,865217, 1004283,947550,1123860,1004274,1004264,956208,956206,1004279,1004296,1004294,947147,1125279,947530,947529,943676,1125668,901914,901918,836537,1124572,836376,905489,972260,972259,968124,1124175,1124819,861006,861008,1124202,1011537,861032,954774,954782,954855,954855,1011061,998485, 1003988,1008663,919080, 1004000, 657749, 931675, 1564232, 1446784,1446799,1446782,581041,581040,1292359,1292360,1292358,585064, 581355, 581356, 579300, 579301,1446821, 1446814,1446808,1317370,1655784, 1651442, 1651210, 1651279, 1217421, 1124489, 1952089,1933067,1933040,1933008,489110,1921689,1937776,1951739,1937835,1937782,1921696,1922017,1922018,1286433,406643,405638, 1974678,2033112,1765829,1771790,1773577,534422,525083,523767,536828,1124344,1124885,1234578,1217445,1222124,1214332,1231553,1240552,1217899,1216096,1214740,1652505,1653657,1653928,1654752,1654964,1652067,1651920,1651563,1651562,1651588,1651561,1651560,1650905,1651636,1651622,1650900,1651623,1654287,1653114,1651758,1124708,1124858,596818,596737,596776,1339767,1125490,1317024,1317024,2025559,1341032,1360578,1118307,1125428,1124426,2335993,1127698,1130925,1141162,1153547,1174787,1169726,1150801,1141097,1151369,1151436,1141094,1151414,1158969,1174147,1166467,1153811,1158967,1166698,1153768,1166858,1151486,1153669,1153671,1133366,1162479,1153536,1151469,1151478,1174067,1151402,2335842,1298701,1300305,581429,581297,583791,580668,580667,2020452,797388,1124701,1450641,1125043,1124683,1444884,578300,578377,578389,580495,2027287,1125087,1125090)) %>% 
  left_join(EVA_ReSu,'plot_id') %>% 
  select(lon, lat) %>% unique() %>% 
  left_join(EVA_ReSu %>% select(plot_id, lon, lat), c('lon','lat')) %>%
  drop_na() %>%
  st_as_sf(coords=c('lon','lat'), crs='WGS84', remove = F)
#mapview::mapview(additional2remove)

EVA_ReSu_spat <- EVA_ReSu_spat %>% 
  anti_join(out_Ukraine %>% st_drop_geometry(), 'plot_id') %>% 
  anti_join(out_Turkey %>% st_drop_geometry(), 'plot_id') %>% 
  anti_join(out_Bulgaria %>% st_drop_geometry(), 'plot_id') %>%
  anti_join(additional2remove %>% st_drop_geometry(), 'plot_id')

# #check again 
st = Sys.time()
out_of_bounds_1km <- st_disjoint(EVA_ReSu_spat, EU_1, sparse = F)
out_of_bounds_1km <- EVA_ReSu_spat[out_of_bounds_1km[,1], ]
Sys.time()-st
#mapview::mapview(out_of_bounds_1km)

# 5. Extract elevation ####

## Load elevation raster data
ele <- paste0('C:/Users/',username,'/OneDrive - CZU v Praze/brno/brno_postdoc/climdata/topodata/elevation/Copernicus_GLO.90_DEM_Europe.tif') %>%
  rast()

## add elevation to plot (take a 2-3 minutes)
st=Sys.time()
EVA_ReSu_spat$elev <- extract(ele, st_transform(EVA_ReSu_spat, crs = crs(ele, proj=T)))$r1
Sys.time()-st

EVA_ReSu_spat %>% # example plot of elevation distribution
  sample_n(5000)%>%
  ggplot(aes(col=elev)) +
  geom_sf() +
  scale_color_viridis_c()


# 6. Assign to either EVA or ReSurveyEU ####

## Transform it back to EVA_ReSu
EVA_ReSu <- EVA_ReSu_spat %>%
  as_Spatial() %>%
  as.data.frame() %>%
  rename(x = coords.x1 , y =  coords.x2)

## Let's re-arrange a bit the order of the variables
EVA_ReSu <- EVA_ReSu %>%
  select(plot_id, dataset, country, lon, lat, x, y, habitat, year,everything())

## Now, subset data from ReSurveyEU...
ResEU <- EVA_ReSu %>% 
  filter(ReSurvey_y_n=='Y') %>%
  filter(manipulate_y_n == 'N') %>%
  # name database
  mutate(database = 'ReSurveyEU', .before = plot_id)

## ... and subset data from EVA...
EVA <- EVA_ReSu %>% 
  filter(ReSurvey_y_n != 'Y' | is.na(ReSurvey_y_n)) %>%
  filter(for_EVA_y_n == 'Y' | is.na(for_EVA_y_n)) %>%
  # name database
  mutate(database = 'EVA', .before = plot_id)

# percentages of plots loss from splitting EVA and ReSurveyEU
(1 - ((nrow(EVA) + nrow(ResEU))/ nrow(EVA_ReSu)))*100

## ... and put them back together.
EVA_ReSu <- bind_rows(EVA, ResEU)
head(EVA_ReSu)
nrow(EVA_ReSu)
table(EVA_ReSu$database)

ggplot() + 
  geom_sf(data = EU, fill='#e6e2df') +
  geom_hex(data=EVA_ReSu, aes(x, y), bins = 70, alpha = 0.8) +
  scale_fill_viridis_c(trans='log10', option='C') +
  facet_wrap(~database) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  labs(fill = 'N. of plots')


# 7. Prepare species data ####
 
## Load core EVA species data
species_data <- read_delim(
  './data/eva/222_PlantDiversity20241017_notJUICE/222_PlantDiversity20241017_notJUICE_species.csv',
  col_select = c(1,3,6,9,10), show_col_types = F
) %>%
  bind_rows( ## add GVRD species data
    read_delim(
      './data/eva/222_PlantDIversity20241025_GVRD/222_GVRD_20241025_notJUICE_species.csv',
      col_select = c(1,3,6,9,10), show_col_types = F
    ) 
  ) %>%
  bind_rows( ## add UKFloodplainMeadows species data
    read_delim(
      './data/eva/222_UKFloodplainMeadows20241031/222_UKFloodplainMeadows20241031_notJUICE_species.csv',
      col_select = c(1,3,6,9,10), show_col_types = F
    ) 
  ) %>%
  setNames(c('plot_id','taxongroup', 'matched_concept', 'layer', 'cover')) %>% # rename columns
  semi_join(EVA_ReSu, 'plot_id') # filter out only the plots we selected this far

# add artificial cover (1) to denmark Nature data
DKnature_plots <- EVA_ReSu %>% filter(dataset=='Denmark Naturdata') %>% select(plot_id)
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

# apply fischer formula to cover data (and remove plots with NA in cover)
fischer.formula <- function(x) {
  1 - (exp(1)^sum(log(1 - x)))
}
st = Sys.time()
species_data_fischer_cover <- species_data %>%
  anti_join(
    DKnature_plots
  ) %>%
  group_by(plot_id) %>% mutate(total.cover = sum(cover)) %>% ungroup() %>%
  group_by(plot_id, matched_concept) %>%
  summarise(rel.cover = sum(cover) / total.cover) %>% # sum the cover of species names
  summarise(rel.cover.fischer = fischer.formula(rel.cover), n.matched = n()) %>% # apply Fischer formula
  ungroup
print(Sys.time() - st)

# Load EIVE
eive <- readxl::read_excel(paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code/data/eiv.data.nomenclature/EIVE.xlsx'), sheet=2) %>%
  select(-contains('.n')) %>% rename(matched_concept=TaxonConcept) 
names(eive)[str_detect(names(eive),'EIVE')] <- paste0('EIV_', names(eive)[str_detect(names(eive),'EIVE')] %>% str_remove_all('EIVEres-'))
eive <- eive %>% select(matched_concept,UUID, EIV_M:EIV_T)

# species list present in EVA/ReSurveyEU (target)
species_list <- species_data %>% 
  select(matched_concept) %>%
  filter(!str_detect(matched_concept, ' species')) %>% # remove genus-only species, we cannot use them for bioindication anyway
  group_by_all() %>%
  summarize(freq=n())
# species_list$new_matched_concept <- species_list$matched_concept
# species_list$new_matched_concept <- ifelse(species_list$new_matched_concept=='Lactuca muralis','Cicerbita muralis',species_list$new_matched_concept)
# species_list$new_matched_concept <- ifelse(species_list$new_matched_concept=='Cerastium fontanum subsp. holosteoides','Cerastium fontanum',species_list$new_matched_concept)
# species_list$new_matched_concept <- ifelse(species_list$new_matched_concept=='Myosotis scorpioides aggr.','Myosotis scorpioides',species_list$new_matched_concept)
# species_list$new_matched_concept <- ifelse(species_list$new_matched_concept=='Bidens tripartitus','Bidens tripartita',species_list$new_matched_concept)
# species_list$new_matched_concept <- ifelse(species_list$new_matched_concept=='Ranunculus acris aggr.','Ranunculus acris',species_list$new_matched_concept)

species_list_j <- species_list %>%
  left_join(eive, 'matched_concept')

# Select species that were not merged
failed_join <- species_list_j %>%
  filter(is.na(UUID)) %>%
  arrange(desc(freq))

#write_csv(failed_join %>% select(1:2) %>% rename(No.EVA.plots=freq), paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code/data/eiv.data.nomenclature/missingEVA.in.EIVE.join.csv'))

##load irena match
irena_main <- readxl::read_excel(
  paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code/data/eiv.data.nomenclature/EIVE_2025-04-16-matched-with-EVA.xlsx'),
  sheet = 1
)
irena_mtch <- readxl::read_excel(
  paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code/data/eiv.data.nomenclature/EIVE_2025-04-16-matched-with-EVA.xlsx'),
  sheet = 2
) %>%
  semi_join(failed_join)
irena_mtch$EIVE %>% table() # 2878 taxa not matched at all with eive
irena_mtch <- irena_mtch %>%
  filter(EIVE=='ok') %>%
  mutate(samename = NameToMatch == matched_concept)

eive_direct <-irena_mtch %>%
  left_join(
    irena_main %>% select(NameToMatch, TaxonConcept), 'NameToMatch'
  ) %>%
  left_join(
    eive %>% select(-UUID) %>% rename(TaxonConcept=matched_concept)
  ) %>%
  arrange(matched_concept) %>%
  split(.,.$matched_concept)

# synthesise EIVs into new values for taxon concept that were not merged:
eive_aggr <- list()
for (i in names(eive_direct)) {
  if(str_detect(i, ' aggr.')){
    
    if(str_detect(i, ' aggr. cf.')){
      if(str_remove(i, ' cf.') %in% eive_direct[[i]]$TaxonConcept){
        dout <- eive_direct[[i]][which(str_remove(i, ' cf.') == eive_direct[[i]]$TaxonConcept),] %>%
          select(matched_concept, contains('EIV_'))
      } else {
        dout <- eive_direct[[i]] %>%
          select(matched_concept, contains('EIV_')) %>%
          group_by(matched_concept) %>%
          summarise_if(is.numeric, median, na.rm=T)
      }

    } else {
      if(str_remove(i, ' aggr.') %in% eive_direct[[i]]$TaxonConcept){
        dout <- eive_direct[[i]][which(str_remove(i, ' aggr.') == eive_direct[[i]]$TaxonConcept),] %>%
          select(matched_concept, contains('EIV_'))
      } else {
        dout <- eive_direct[[i]] %>%
          select(matched_concept, contains('EIV_')) %>%
          group_by(matched_concept) %>%
          summarise_if(is.numeric, median, na.rm=T)
      }
    }
  
  } else {
    dout <- eive_direct[[i]] %>%
      select(matched_concept, contains('EIV_')) %>%
      group_by(matched_concept) %>%
      summarise_if(is.numeric, median, na.rm=T) 
  }
  eive_aggr[[i]] <- dout
}
eive_aggr <- eive_aggr %>% bind_rows()

# new join!
eive_merge <- eive %>%
  select(-UUID) %>%
  bind_rows(eive_aggr)
species_list_nj <- species_list %>%
  left_join(eive_merge, 'matched_concept')

# will the species be included in the analysis?
species_list_nj <- species_list_nj %>%
  mutate(EIVE.match = !if_all(contains('EIV_'), is.na), .after = matched_concept)
table(species_list_nj$EIVE.match)
# FALSE  TRUE 
#  2936 13874 
round(prop.table(table(species_list_nj$EIVE.match)),2)
# FALSE  TRUE 
# 0.17  0.83 

failed_taxa_final <- species_list_nj %>% filter(!EIVE.match)
failed_taxa_final %>% pull(freq) %>% range() # from 1 to 1014
failed_taxa_final %>% pull(freq) %>% mean() # 18.2 species
failed_taxa_final %>% pull(freq) %>% sd() #53.6 species

EIVE_EVA_merge_filename <- paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code/data/eiv.data.nomenclature/EIVE_EVA_merge_', date.update, '.csv')
write_csv(species_list_nj %>% filter(EIVE.match) %>% select(-EIVE.match), 
          file = EIVE_EVA_merge_filename)

# 8. Calculate CMs ####

species_data_eive <-  species_data %>%
  {if(herbsOnly) anti_join(., drevojan_treeshrub, by='matched_concept') else .} %>%
  select(plot_id, matched_concept) %>%
  ungroup() %>%
  unique() %>%
  left_join(eive_merge)

glimpse(species_data_eive)

# calc cm (community means)
n_na <- \(x){sum(rep(1, length(x))[!is.na(x)])}
mean_na <- \(x){mean(x, na.rm = T)}
sd_na <- \(x){sd(x, na.rm = T)}

st=Sys.time() # approx. 1.5 mins
CM_eive <- species_data_eive %>%
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
hist(CM_eive$cm.EIV_M)
hist(CM_eive$cm.EIV_N)
hist(CM_eive$cm.EIV_R)
hist(CM_eive$cm.EIV_L)
hist(CM_eive$cm.EIV_T)

table(CM_eive$n.EIV_M >= 0.8) %>% prop.table()
table(CM_eive$n.EIV_N >= 0.8) %>% prop.table()
table(CM_eive$n.EIV_R >= 0.8) %>% prop.table()
table(CM_eive$n.EIV_L >= 0.8) %>% prop.table()
table(CM_eive$n.EIV_T >= 0.8) %>% prop.table()

table(CM_eive$plot_id %in% EVA_ReSu$plot_id)
table(EVA_ReSu$plot_id %in% CM_eive$plot_id ) # there are a few plots without any vascular plants (510)
# plots_with_NA_CM <- EVA_ReSu[!(EVA_ReSu$plot_id %in% CM_eive$plot_id),]

## Only select plots with actual calculable CMs
EVA_ReSu <- EVA_ReSu[(EVA_ReSu$plot_id %in% CM_eive$plot_id),]

## discard all plots where none of the indicator CM has relative n >= 0.8
CM_eive_low_treshold_all = CM_eive[!(CM_eive$n.EIV_M >= 0.8 | CM_eive$n.EIV_N >= 0.8 | CM_eive$n.EIV_T >= 0.8 | CM_eive$n.EIV_L >= 0.8 | CM_eive$n.EIV_R >= 0.8 ),]
CM_eive <- CM_eive %>% anti_join(CM_eive_low_treshold_all, 'plot_id')
EVA_ReSu <- EVA_ReSu %>% anti_join(CM_eive_low_treshold_all, 'plot_id')

table(EVA_ReSu$plot_id %in% CM_eive$plot_id )

## Load richness outliers, checked by I. Knollova and M. Chytry
richness_outliers <- 
  paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU/data/high.richnes.outliers_IlonaKnollova.xlsx') %>%
  readxl::read_excel() %>%
  select(1) %>%
  setNames(c('plot_id'))

## Remove richness outliers
CM_eive <- CM_eive %>%
  anti_join(
    richness_outliers
  )

EVA_ReSu <- EVA_ReSu %>% 
  semi_join(CM_eive, 'plot_id') %>% # filter header data accordingly
  left_join(CM_eive) # add data

## remove plots with species richness higher than 110
dat_high_S_check <- EVA_ReSu %>%
  filter(n >= 110)
EVA_ReSu <- EVA_ReSu %>%
  anti_join(dat_high_S_check, 'plot_id')


#9. Finalize main data & export ####

EVA_ReSu <- EVA_ReSu %>%
  filter(country != 'Turkey') # remove unwanted countries, ca 7 plots are left and labelled as Turkey

# filter plots for which no species with EIVE values are found
EVA_ReSu_all.missing.EIVE <- EVA_ReSu %>%
  filter((is.na(cm.EIV_M) & is.na(cm.EIV_N) & is.na(cm.EIV_R) & is.na(cm.EIV_L) & is.na(cm.EIV_T))) 
EVA_ReSu_all.missing.EIVE %>%
  nrow() # currently 66 in total

EVA_ReSu <- EVA_ReSu %>% 
  anti_join(EVA_ReSu_all.missing.EIVE, 'plot_id')

## Finalize EVA data
EVA <- EVA_ReSu %>%
  filter(database == 'EVA') %>%
  select(
    database, plot_id, dataset, ESy, habitat, lon, lat, x, y, elev, year, plot_size 
  ) %>%
  arrange(plot_id) %>%
  left_join(CM_eive, 'plot_id')

## Finalize ReSurveyEU data
ResEU <- EVA_ReSu %>%
  filter(database == 'ReSurveyEU') %>%
  select(
    database, plot_id, dataset, ESy, habitat, lon, lat, x, y, elev, year, plot_size,
    ReSurvey_type, ReSur_site, ReSur_plot, ReSur_obs, ReSur_time
  ) %>%
  left_join(CM_eive, 'plot_id')
ResEU$ReSurvey_type %>% table()
ResEU$ReSur_type <- ifelse(str_detect(ResEU$ReSurvey_type, 'esampling'), 'Resampling', 'Permanent')
ResEU$ReSur_type %>% table() %>% prop.table()
ResEU <- ResEU %>%
  select(database, plot_id, contains('ReSur'), everything()) %>%
  arrange(plot_id) %>%
  select(-ReSurvey_type)

## Export EVA data
write_csv(EVA, paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code', '/data/input/EVA', '_', date.update,'.csv'))
## Export ReSurveyEU data
write_csv(ResEU, paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code', '/data/input/ReSurveyEU','_', date.update,'.csv'))


# Pretty maps
library(tidyverse);library(sf);library(terra)
username <- str_split(getwd(), '/')[[1]][3]
setwd(paste0('C:/Users/', username, '/OneDrive - CZU v Praze/czu/intrplEU_EIV/'))

prj = 25832 
EU <- read_sf(list.files(paste0('C:/Users/',username,'/OneDrive - CZU v Praze/brno/brno_postdoc/maps/euro+med_map/Final.Map'), pattern = 'shp',full.names = T)) %>%
  semi_join(data.frame(name =
                         c('Albania', 'Austria', 'Baleares', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria',
                           'Corsica', 'Crete', 'Croatia', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany',
                           'Greece', 'Hungary', 'Ireland', 'Italy', 'Kaliningrad Region', 'Kosovo', 'Latvia', 'Liechtenstein',
                           'Lithuania', 'Luxembourg', 'Malta', 'Moldova', 'Montenegro', 'Netherlands', 'North Macedonia',
                           'Norway', 'Poland', 'Portugal', 'Romania', 'Sardinia', 'Serbia', 'Sicily',
                           'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine', 'United Kingdom'
                         )), 'name') %>%
  st_transform(crs = prj) %>%
  st_simplify(dTolerance = 1000)
EU <- EU %>% st_buffer(1*1000)
EU <- st_as_sf(st_union(EU$geometry)) 

d2p <- bind_rows(ResEU %>% select(database, x, y), EVA %>% select(database, x, y)) 

p <- ggplot() + 
  geom_sf(data = EU, fill='#e6e2df') +
  geom_hex(data=d2p, aes(x, y), bins = 70, alpha = 0.8) +
  scale_fill_viridis_c(trans='log10', option='C') +
  facet_wrap(~database) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  labs(fill = 'N. of plots')
p

#ggsave('./fig/lon_lat_plot.pdf', p, width = 8, height = 4.6)
ggsave('./fig/lon_lat_plot.jpg', p, width = 8, height = 4.6, dpi = 600)

grd <- st_make_grid(st_bbox(EU), cellsize=50*1000) %>%
  st_as_sf() %>%
  mutate(grd.id=1:nrow(.))

d2p <- bind_rows(ResEU %>% select(plot_id, database, x, y, contains('cm.EIV')), EVA %>% select(plot_id, database, x, y, contains('cm.EIV')))%>%
  st_as_sf(
    coords=c('x','y'), crs=crs(EU)
  )

grd <- st_join(grd,d2p)

mean_na <- \(x){mean(x, na.rm = T)}

st = Sys.time()
grd <- grd %>% 
  group_by(grd.id) %>%
  select(-plot_id, -database) %>%
  summarise_if(
    is.numeric,
    mean_na
  )
Sys.time() - st

grd_tidy <- grd %>%
  ungroup() %>%
  gather('var','val',contains('cm.EIV'))

for (i in c('cm.EIV_M', 'cm.EIV_N', 'cm.EIV_R', 'cm.EIV_L', 'cm.EIV_T')) {
  d2p <- grd_tidy %>% filter(var==i) %>% filter(!is.na(val))
  p_CWM <- ggplot() +
    geom_sf(data = EU, fill='#e6e2df') +
    geom_sf(data = d2p, aes(fill=val), color=NA, alpha = 0.8) +
    scale_fill_viridis_c(option='C') +
    theme_bw() +
    theme(axis.title = element_blank()) +
    labs(fill = str_remove(i, 'cm.')) + ggtitle(str_remove(i, 'cm.'))
  p_CWM
  
  ggsave(paste0('./fig/50km.summary_',str_remove(i, 'cm.'),'.pdf'), p_CWM, width = 4.2, height = 4.6)
  
}

#10. Prepare supplementary data for community weighted means ####

## Calc CWM (weighted mean for supplementary analysis)
w_mean_na <- \(x_val, w_val){weighted.mean(x_val, w_val, na.rm = T)}
st = Sys.time()
CWM_eive <- species_data_fischer_cover %>% select(-n.matched) %>%
  semi_join(bind_rows(EVA, ResEU), 'plot_id') %>%
  left_join(eive_merge, 'matched_concept') %>%
  group_by(plot_id) %>%
  summarise(
    cwm.EIV_M=w_mean_na(EIV_M, w_val = rel.cover.fischer),
    cwm.EIV_N=w_mean_na(EIV_N, w_val = rel.cover.fischer),
    cwm.EIV_R=w_mean_na(EIV_R, w_val = rel.cover.fischer),
    cwm.EIV_L=w_mean_na(EIV_L, w_val = rel.cover.fischer),
    cwm.EIV_T=w_mean_na(EIV_T, w_val = rel.cover.fischer)
  ) %>%
  ungroup()
Sys.time()-st

EVA_ReSu_CWM <- CWM_eive %>%
  semi_join(EVA_ReSu) %>%
  left_join(EVA_ReSu) %>%
  select(plot_id, contains('m.EIV_'), everything())

## export
write_csv(EVA_ReSu_CWM, paste0('C:/Users/',username,'/OneDrive - CZU v Praze/czu/intrplEU_EIV/code', '/data/input/EVA_ReSu_CWM', '_', date.update,'.csv'))