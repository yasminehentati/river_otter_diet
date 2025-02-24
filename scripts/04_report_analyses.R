
# packages
p_load(corrplot,glmmTMB,AICcmodavg,sjPlot,DHARMa,MuMIn,ggplot2,janitor,dplyr,tidyr,
       readr,here,viridis)


####################################
##### first we'll add some additional latrine info 
# let's classify latrines by freshwater/saltwater 

### aquatic 
# read in water
water <- st_read(here("spatial_data", 
                      "DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp")) 

sf_use_s2(FALSE)
st_crs(water)

# crop to our map area and put in the projection we're working in
water <-  water %>%
  st_transform(crs = 4326) %>% 
  st_crop(c(xmin= -121.7, ymin = 46.7, xmax = -122.8, ymax = 48.5)) %>%
  st_transform(crs = 32610) 

# read in our site locations 
sites <- st_read(here("clean_data", "latrines.shp")) %>% 
  st_transform(crs = "EPSG:32610")



# let's add a new attribute to water that will be based on hydrography categories 
# according to the metadata https://fortress.wa.gov/dnr/adminsa/GisData/metadata/ROPA_WBHYDRO.pdf
# and the shoerline data set https://geo.wa.gov/datasets/wadnr::shorezone-inventory-salt-marsh/explore?location=47.903540%2C-122.176510%2C10.00
# all saltwater bodies are encompassed under "bay"

water <- water %>%
  mutate(water_type = ifelse(WB_HYDR_FT == "Bay", "saltwater", "freshwater"))


# now let's calculate distance to saltwater for each site
# this is a straight line closest distance 

saltwater <- water %>% filter(water_type == "saltwater")

# calculate the nearest distance from each site to the closest saltwater feature
sites <- sites %>%
  mutate(saltwater_dist = st_distance(geometry, st_geometry(saltwater)) %>%
           apply(1, min))  # get the minimum distance for each point

print(sites)




# now let's concatenate sites with our other spatial covariates 
 

all_covs <- read_csv(here("spatial_data", "covariates", "ALL_ENV_URB_SITES_2000m.csv"))


site_covs <- left_join(sites, all_covs, by = c("city", "site")) %>%
select(-c(1:5,8:13))

colnames(site_covs)



write_csv(site_covs, here("spatial_data", "covariates", "ALL_ENV_URB_WA_SITES_2000m.csv"))


############## 


##############
### let's create  diet classification columns 
# first we'll read in the diet data 

verts <- read_csv(here("clean_data", "12S_final.csv"))
inverts <- read_csv(here("clean_data", "16S_final.csv"))
site_covs <- read_csv(here("spatial_data", "covariates", "ALL_ENV_URB_WA_SITES_2000m.csv"))



#################
# diet classifications 
# for freshwater/saltwater, checked using https://www.sealifebase.se/
## freshwater malacostracae
# Asellidae (isopod), Cambaridae (crayfish), Astacidae (crayfish)

fresh_mala_fam <- c("Asellidae", "Cambaridae", "Astacidae")

inverts <- inverts %>%
  mutate(diet_class = ifelse(family %in% fresh_mala_fam, "fresh_mala", NA))

sum(inverts$diet_class == "fresh_mala", na.rm = TRUE)



## saltwater malacostracae
# MANUAL LIST: 
# Amphiphoidae (amphipods), Aoridae (amphipod), Cancridae (crabs), Caprellidae (amphipod), 
# Corophiidae (amphipod -- salt and fresh), Dulichiidae (amphipod), Ischyroceridae (amphipod), 
# Sphaeromatidae (isopod), Varunidae (crabs), Panopeidae (mud crabs), Oregoniidae (crabs), 
# Paguridae (hermit crabs), Pandalidae (shrimps), Cheiragonidae (helmet crab)


inverts %>%
  filter(class == "Malacostraca") %>%
  pull(family) %>%
  unique()

salt_mala_fam <- c("Amphiphoidae", "Aoridae", "Cancridae", "Caprellidae", 
                         "Corophiidae", "Dulichiidae", "Ischyroceridae", "Sphaeromatidae", 
                         "Varunidae", "Panopeidae", "Oregoniidae", "Paguridae", 
                         "Pandalidae", "Cheiragonidae")

inverts <- inverts %>%
  mutate(diet_class = ifelse(family %in% salt_mala_fam, "salt_mala", diet_class))

sum(inverts$diet_class == "salt_mala", na.rm = TRUE)

## dungeness crabs (overwrite dungeness in saltwater mala)
#dungeness <- c("Metacarcinus") 

#inverts <- inverts %>%
#  mutate(diet_class = ifelse(genus %in% dungeness, "dungeness", diet_class))

#sum(inverts$diet_class == "dungeness", na.rm = TRUE)

## insects
insects <- c("Insecta") 

inverts <- inverts %>%
  mutate(diet_class = ifelse(class %in% insects, "insects", diet_class))

sum(inverts$diet_class == "insects", na.rm = TRUE)

## gastropods
gastropods <- c("Gastropoda") 

inverts <- inverts %>%
  mutate(diet_class = ifelse(class %in% gastropods, "gastropods", diet_class))

sum(inverts$diet_class == "gastropods", na.rm = TRUE)

#########
# vert classifications 

## ray finned fishes 

# centrarchiformes, perciformes, blenniiformes (perch-like + surf perch)
# cypriniformes (cyprinid fish -- carps and minnows)
# siluriformes (catfish -- all brown bullhead)
# clupeiformes (pacific herring -- not enough to isolate)
# pleuronectiformes (flatfish -- 2 spp)
# salmoniformes 

verts$class

fish <- c("Actinopteri") 

verts <- verts %>%
  mutate(diet_class = ifelse(class %in% fish, "fish", NA))

sum(verts$diet_class == "fish", na.rm = TRUE)

verts %>%
  filter(class == "Actinopteri") %>%
  pull(order) %>%
  unique()

head(all_diet)

# breaking down groups of ray-finned fish

# perch-like 
perch <- c("Blenniiformes", "Perciformes", "Cantrarchiformes") 

verts <- verts %>%
  mutate(diet_class = ifelse(order %in% perch, "perch", diet_class))

sum(verts$diet_class == "perch", na.rm = TRUE)

# cypriniformes
cypri <- c("Cypriniformes") 

verts <- verts %>%
  mutate(diet_class = ifelse(order %in% cypri, "cypri", diet_class))

sum(verts$diet_class == "cypri", na.rm = TRUE)

# siluriformes (brownbullhead)
bullhead <- c("Siluriformes") 

verts <- verts %>%
  mutate(diet_class = ifelse(order %in% bullhead, "bullhead", diet_class))

sum(verts$diet_class == "bullhead", na.rm = TRUE)

# pleuronectiformes
flatfish <- c("Pleuronectiformes") 

verts <- verts %>%
  mutate(diet_class = ifelse(order %in% flatfish, "flatfish", diet_class))

sum(verts$diet_class == "flatfish", na.rm = TRUE)

unique(verts$diet_class)

## ducks  & rails 
verts %>%
  filter(class == "Aves") %>%
  pull(family) %>%
  unique()

ducksrails <- c("Anatidae", "Rallidae") 

verts <- verts %>%
  mutate(diet_class = ifelse(family %in% ducksrails, "ducksrails", diet_class))

sum(verts$diet_class == "ducksrails", na.rm = TRUE)

## passerine birds 

passerine <- c("Corvidae", "Turdidae") 

verts <- verts %>%
  mutate(diet_class = ifelse(family %in% passerine, "passerine", diet_class))

sum(verts$diet_class == "passerine", na.rm = TRUE)

## small mammals -- non-anthro
verts %>%
  filter(class == "Mammalia") %>%
  pull(family) %>%
  unique()

smallmammals <- c("Muridae", "Castoridae", "Leporidae", "Cricetidae")

verts <- verts %>%
  mutate(diet_class = ifelse(family %in% smallmammals, "smallmammals", diet_class))

sum(verts$diet_class == "smallmammals", na.rm = TRUE)

## anthropogenic mammals 
# turkey, chicken, sheep, cow, cat, pig 
anthro <- c("Bos" , "Phasianidae" , "Oves" , "Gallus", "Felis", "Sus")

verts <- verts %>%
  mutate(diet_class = ifelse(genus %in% anthro, "anthro", diet_class))

sum(verts$diet_class == "anthro", na.rm = TRUE)

## salmonids (overwrite salmonids in ray finned)

salmon <- c("Salmonidae") 

verts <- verts %>%
  mutate(diet_class = ifelse(family %in% salmon, "salmonids", diet_class))

sum(verts$diet_class == "salmonids", na.rm = TRUE)


############## 
## combine all data 
inverts$processed_taxon
inverts <- inverts %>% 
  select(-possible_taxon,superclass)

# Get the full set of column names from both data frames
all_columns <- union(colnames(verts), colnames(inverts))

# Ensure verts has all columns by adding missing ones with NA
for (col in setdiff(all_columns, colnames(verts))) {
  verts[[col]] <- NA
}

# Ensure inverts has all columns by adding missing ones with NA
for (col in setdiff(all_columns, colnames(inverts))) {
  inverts[[col]] <- NA
}

# Reorder columns to match the full set
verts <- verts %>% select(all_columns)
inverts <- inverts %>% select(all_columns)

# Bind the two data frames together
all_diet <- as.data.frame(bind_rows(verts, inverts))
class(all_diet)


#write_csv(all_diet, here("clean_data", "combined_diet_16s_12s.csv"))

# let's add covariates to diet data 
# first we need to add site names to the diet list 

## re-read in all the data 

all_diet <- read_csv(here("clean_data", "combined_diet_16s_12s.csv"))
site_covs <- read_csv(here("spatial_data", "covariates", "ALL_ENV_URB_WA_SITES_2000m.csv"))

metadata <- read_csv(here("data_from_ellie", "field_site_data.csv")) %>% 
  rename(sample_ID = "Scat ID",
         site = "Site",
         date = "Date", 
         month = "Month", 
         city = "City") %>% 
  select(sample_ID, site, date, month, city)
  
all_diet <- left_join(all_diet, metadata, by = c("sample_ID"))

# now join the spatial data


all_diet <- left_join(all_diet, site_covs, by = c("city", "site"))


# for some reason LCP047 is missing site data 
# it's from FGL so we'll fill the columns based on that 
colnames(all_diet)

all_diet <- all_diet %>%
  mutate(site = ifelse(grepl("LC047", sample_ID), "FGL", site)) %>%
  group_by(site) %>%
  mutate(across(35:62, ~ ifelse(sample_ID == "LC047" & is.na(.), first(na.omit(.)), .))) %>%
  ungroup()

# finally let's add water body names to each site before we run glmmms


all_diet <- all_diet %>% 
  mutate(waterbody = case_when(
    site %in% c("LNM", "LWM", "MBB", "LSM", "LBP", "WNP", "CBM") ~ "Lake Washington",
    site %in% c("FGL", "DSS", "CWG", "BRF") ~ "Duwamish/Green River",
    site %in% c("BSW", "WPA") ~ "Portage & Union Bays",
    site %in% c("Cornet Bay", "GGP", "DNB") ~ "Puget Sound", 
    site %in% c("RVL", "GGC") ~ "Kelsey Creek")
  )



# let's make our variables factors 
all_diet$waterbody <- as.factor(all_diet$waterbody)
all_diet$site <- as.factor(all_diet$site)

write_csv(all_diet, here("clean_data", "combined_diet_covs_all.csv"))

# check variable correlatoin 
cordat <- all_diet %>% distinct(site, .keep_all = TRUE)  %>% 
  dplyr::select(c(saltwater_dist, POPDEN2020, RISK,50:60))

corr <- cor(cordat, method = "spearman", use="pairwise.complete.obs") 

corrplot(corr, 
         method = 'color', 
         type = 'lower', 
         order = 'original', 
         tl.col = 'black',
         cl.ratio = 0.2, 
         tl.srt = 45, 
         col = viridis(200),  # Use viridis for colorblind-friendly palette
         cl.pos = 'b', 
         addgrid.col = 'white', 
         addCoef.col = 'black')  # Coefficient text color

unique(all_diet$diet_class)

nrow(all_diet)


# quick summary 
diet_summary <- all_diet %>%
  summarise(
    num_unique_samples = n_distinct(sample_ID),                      # Unique sample_IDs
    avg_reads_per_sample = mean(total_reads_avg, na.rm = TRUE),      # Mean of total_reads_avg
    num_unique_families = n_distinct(family, na.rm = TRUE),          # Unique families
    num_unique_species = n_distinct(species, na.rm = TRUE)           # Unique species
  )

print(diet_summary)

####################
# glmmms 

# first we need to organize the data into presence/absence by making new binary columns 
# for each of our diet classes of interest 

all_diet <- all_diet %>%
  mutate(
    fish = ifelse(diet_class == "fish", 1, ifelse(is.na(diet_class), 0, 0)),
    ducksrails = ifelse(diet_class == "ducksrails", 1, ifelse(is.na(diet_class), 0, 0)),
    anthro = ifelse(diet_class == "anthro", 1, ifelse(is.na(diet_class), 0, 0)),
    smallmammals = ifelse(diet_class == "smallmammals", 1, ifelse(is.na(diet_class), 0, 0)),
    salmonids = ifelse(diet_class == "salmonids", 1, ifelse(is.na(diet_class), 0, 0)),
    passerine = ifelse(diet_class == "passerine", 1, ifelse(is.na(diet_class), 0, 0)),
    fresh_mala = ifelse(diet_class == "fresh_mala", 1, ifelse(is.na(diet_class), 0, 0)),
    salt_mala = ifelse(diet_class == "salt_mala", 1, ifelse(is.na(diet_class), 0, 0)),
    dungeness = ifelse(diet_class == "dungeness", 1, ifelse(is.na(diet_class), 0, 0)),
    perch = ifelse(diet_class == "perch", 1, ifelse(is.na(diet_class), 0, 0)),
    cypri = ifelse(diet_class == "cypri", 1, ifelse(is.na(diet_class), 0, 0)),
    bullhead = ifelse(diet_class == "bullhead", 1, ifelse(is.na(diet_class), 0, 0)),
    flatfish = ifelse(diet_class == "flatfish", 1, ifelse(is.na(diet_class), 0, 0)),
    insects = ifelse(diet_class == "insects", 1, ifelse(is.na(diet_class), 0, 0)),
    gastropods = ifelse(diet_class == "gastropods", 1, ifelse(is.na(diet_class), 0, 0)),
  ) %>% 
  mutate(across(starts_with("fish"):starts_with("gastropods"), ~ replace(., is.na(.), 0)))  

View(all_diet)
# let's save our final data set 
write_csv(all_diet, here("clean_data", "all_diet_cov_data.csv"))

##### 
# freshwater mala 
fresh_mala <- glmmTMB(fresh_mala ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                        (1|waterbody/site), data=all_diet, family="binomial")
summary(fresh_mala)

plot_model(fresh_mala, type="est", title="Freshwater Malacostracans", show.p=T)


# saltwater mala 
salt_mala <- glmmTMB(salt_mala ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                       (1|waterbody/site), data=all_diet, family="binomial")

summary(salt_mala)

# dungeness 
#dungeness <- glmmTMB(dungeness ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
 #                      (1|site), data=all_diet, family="binomial")

#summary(dungeness)

# insects
all_diet$insects
insects <- glmmTMB(insects ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                       (1|waterbody/site), data=all_diet, family="binomial")

summary(insects)

# gastropods
gastropods <- glmmTMB(gastropods ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                     (1|site), data=all_diet, family="binomial")

summary(gastropods)

# ducks/rails 
ducksrails <- glmmTMB(ducksrails ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                       (1|site), data=all_diet, family="binomial")

summary(ducksrails)

# ray finned fishes 
fish <- glmmTMB(fish ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                        (1|site), data=all_diet, family="binomial")
summary(fish)


# perch
perch <- glmmTMB(perch ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                  (1|site), data=all_diet, family="binomial")
summary(perch)

# cyprinid fish
cypri <- glmmTMB(cypri ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                   (1|site), data=all_diet, family="binomial")
summary(cypri)

# bullhead
bullhead <- glmmTMB(bullhead ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                   (1|site), data=all_diet, family="binomial")
summary(bullhead)

# flatfish
flatfish <- glmmTMB(flatfish ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                      (1|site), data=all_diet, family="binomial")
summary(flatfish)

# salmonids 
salmonids <- glmmTMB(salmonids ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                  (1|site), data=all_diet, family="binomial")
summary(salmonids)

# small mammals 
smallmammals <- glmmTMB(smallmammals ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                       (1|site), data=all_diet, family="binomial")
summary(smallmammals)

# anthro 
anthro <- glmmTMB(anthro ~ scale(RISK) + scale(POPDEN2020) + scale(saltwater_dist) + 
                          (1|site), data=all_diet, family="binomial")
summary(anthro)


################ 
# PLOT 

# inverts 
(plot_global <- plot_models(fresh_mala, salt_mala, insects, 
                            show.values = T,
                            show.legend = F,
                            grid = TRUE, vline.color = "gray",
                           #   grid.breaks = c(.1,1,10),
                            p.shape     = "FALSE"#,
                            #    rm.terms = c("") , 
                            #       transform = NULL
                             #  axis.lim = c(0.1,10) 
)+ 
  scale_y_log10(limits = c(-0.1,10)) + 
  scale_color_sjplot("eight") + 
  #  scale_y_continuous(labels = scales::comma_format()) +
  scale_x_discrete(labels=c("Contaminant risk", "Human pop. density", "Dist. to saltwater"), limits = rev) + 
  font_size(axis_title.x = 0, axis_title.y = 12, labels.x = 8, labels.y = 12))

# birds and mammals 
(plot_global <- plot_models(ducksrails,smallmammals, # anthro,
                            show.values = T,
                            show.legend = F,
                            grid = TRUE, vline.color = "gray",
                               #  grid.breaks = c(.1,1,10),
                            p.shape     = "FALSE"#,
                            #    rm.terms = c("") , 
                            #       transform = NULL
                               #    axis.lim = c(0.1,10) 
)+ 
       scale_y_log10(limits = c(0.1,10)) + 
    scale_color_sjplot("eight") + 
  #  scale_y_continuous(labels = scales::comma_format()) +
   scale_x_discrete(labels=c("Contaminant risk", "Human pop. density", "Dist. to saltwater"), limits = rev) + 
    font_size(axis_title.x = 0, axis_title.y = 12, labels.x = 8, labels.y = 12))

# fish 
(plot_global <- plot_models(salmonids,bullhead,cypri,perch,flatfish,
                            show.values = T,
                            show.legend = F,
                            grid = TRUE, vline.color = "gray",
                            #     grid.breaks = c(.1,1,10),
                            p.shape     = "FALSE"#,
                            #    rm.terms = c("") , 
                            #       transform = NULL
                            #       axis.lim = c(0.1,10) 
)+ 
        scale_y_log10(limits = c(-0.1, 10)) + 
    scale_color_sjplot("eight") + 
    scale_x_discrete(labels=c("Contaminant risk", "Human pop. density", "Dist. to saltwater"), limits = rev) + 
    font_size(axis_title.x = 0, axis_title.y = 12, labels.x = 8, labels.y = 12))


#################
# seasonal analysis
# Create new season column based on months
all_diet <- all_diet %>%
  mutate(month,  # Ensure months is numeric
    season = case_when(
      month %in% c(1, 2) ~ "winter",   # December, January, February
      month %in% c(3, 4, 5) ~ "spring",    # March, April, May
      month %in% c(6, 7, 8) ~ "summer",    # June, July, August
      month %in% c(9) ~ "fall",    # September, October, November
      TRUE ~ NA_character_                 # For any other value (if exists)
    )
  )

# now let's keep only spring and summer data set

all_diet_ssn <- all_diet %>%
  filter(season %in% c("summer", "spring"))

# convert to factor 
all_diet_ssn$season <- as.factor(all_diet_ssn$season)

## fit logistic regressions 

fresh_mala_ssn <- glmmTMB(fresh_mala ~ season, data=all_diet_ssn, family="binomial")
summary(fresh_mala_ssn)

plot_model(fresh_mala_ssn, type="est", title="Freshwater Malacostracans", show.values=T)


# saltwater mala 
salt_mala_ssn <- glmmTMB(salt_mala ~ season, data=all_diet_ssn, family="binomial")

summary(salt_mala_ssn)

plot_model(salt_mala_ssn, type="est", title="Saltwater Malacostracans", show.p=T)


# dungeness 
dungeness_ssn <- glmmTMB(dungeness ~ season, data=all_diet_ssn, family="binomial")

summary(dungeness_ssn)

# gastropods
gastropods_ssn <- glmmTMB(gastropods ~ season, data=all_diet_ssn, family="binomial")

summary(gastropods_ssn)

# insects
insects_ssn <- glmmTMB(insects ~ season, data=all_diet_ssn, family="binomial")

summary(insects_ssn)

# ducks/rails 
ducksrails_ssn <- glmmTMB(ducksrails ~ season, data=all_diet_ssn, family="binomial")

summary(ducksrails_ssn)

# ray finned fishes 
fish_ssn <- glmmTMB(fish ~ season, data=all_diet_ssn, family="binomial")

summary(fish_ssn)

# perch
perch_ssn <- glmmTMB(perch ~ season, data=all_diet_ssn, family="binomial")
summary(perch_ssn)

# cyprinid fish
cypri_ssn <- glmmTMB(cypri ~ season, data=all_diet_ssn, family="binomial")
summary(cypri_ssn)

# bullhead
bullhead_ssn <- glmmTMB(bullhead ~ season, data=all_diet_ssn, family="binomial")
summary(bullhead_ssn)

# flatfish
flatfish_ssn <- glmmTMB(flatfish ~ season, data=all_diet_ssn, family="binomial")
summary(flatfish_ssn)

# salmonids 
salmonids_ssn <- glmmTMB(salmonids ~ season, data=all_diet_ssn, family="binomial")
summary(salmonids_ssn)

# small mammals 
smallmammals_ssn <- glmmTMB(smallmammals ~ season, data=all_diet_ssn, family="binomial")
summary(smallmammals_ssn)

# anthro 
anthro_ssn <- glmmTMB(anthro ~ season, data=all_diet_ssn, family="binomial")
summary(anthro_ssn)



tab_model(insects_ssn, fresh_mala_ssn, salt_mala_ssn, 
          collapse.ci = TRUE, 
          p.style     = "numeric_stars")
) 

tab_model(flatfish_ssn, perch_ssn, cypri_ssn, bullhead_ssn, salmonids_ssn,
          collapse.ci = TRUE, 
          p.style     = "numeric_stars")
) 

tab_model(smallmammals_ssn, ducksrails_ssn, anthro_ssn,
          collapse.ci = TRUE, 
          p.style     = "numeric_stars")
) 