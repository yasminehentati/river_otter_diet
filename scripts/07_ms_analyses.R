### MS analyses 


# packages
library(pacman)
p_load(corrplot,glmmTMB,AICcmodavg,sjPlot,DHARMa,MuMIn,ggplot2,janitor,dplyr,tidyr,
       readr,here,viridis)

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



Let's look at potential associations between dietary groups. 
```{r}

## potential associations 
# beetles families and tachnidae 

# birds and bird mites (Arachnida)

# cyprinid fish and ligula intestinalis 


```

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