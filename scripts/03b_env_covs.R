#### spatial covs 2: env health data 

###########################
# income data 
# let's read in some tidycensus data -- this will help us to assign GEOID values for 
# env health data 


## median household income data from American Community Survey
Sys.getenv("CENSUS_API_KEY")
# load census API key
census_api_key("e649b78a1d98fe7e1c9ff7039840781976777eb6",
               install = TRUE)
readRenviron("~/.Renviron")

# load in ACS variables  from 2019 
v19 <- load_variables(2019, "acs5", cache = TRUE)


# load in shapefile of subdivisions  
# with median household income as variable of interest

# for washington 
tractincomeWA <- get_acs(state = "WA", 
                         geography = "tract", 
                         variables = c(medincome="B19013_001"),
                         year = 2019, geometry = TRUE) # need to specify 2015-2019 data 
# since acs now automatically calls 2017-2021 and GEOIDs change each census --
# the env. health data has 2019 GEOIDs 

# cut out everything except the counties we want
# do this using the GEOID codes for each county 

tractsKP <- tractincomeWA %>% dplyr::filter(substr(GEOID, 1, 5) 
                                            %in% c("53033", "53053", # king pierce 
                                                   "53061", "53035", # snoho kitsap
                                                   "53029", "53057")) # island skagit



mapview(tractsKP) # looks good

# crop to actual study area 
#tractsKP <- st_transform(tractsKP, crs=4326)
#tractsKP_crop <- st_crop(tractsKP, c(xmin= -121.7, ymin = 46.7, xmax = -122.8, ymax = 47.8))

# write shapefiles for cropped area 
st_write(tractsKP, here("spatial_data", "w_wa_urban_med_income.shp"),
         append = FALSE)



######### environmental metrics 



# read in env health data 

lead <- read_csv(here("spatial_data", "env_health_data", "Lead_Risk_from_Housing.csv"), 
                 name_repair = janitor::make_clean_names) %>% 
  dplyr::select(census_tract, number_units_w_lead_hazard_estimated)
#        col_select = c("County_Name", "Census_Tract", "Units_Lead", 
#                 "Total_Units", "Pct_Units_Lead", "IBL_Rank"))
colnames(lead)

pm25 <- read_csv(here("spatial_data", "env_health_data", "PM2.5_Concentration.csv"), 
                 name_repair = janitor::make_clean_names) %>% 
  select(census_tract, count)
#    col_names=c("County_Name", "YearGroup", "Census_Tract", "PM2.5_Count", "IBL_Rank"))

prox_haz_waste <- read_csv(here("spatial_data", "env_health_data", 
                                "Proximity_to_Hazardous_Waste_Treatment_Storage_and_Disposal_Facilities_(TSDFs).csv"), 
                           name_repair = janitor::make_clean_names) %>% 
  select(census_tract, average_ptsdf)
# col_names=c("County_Name", "Census_Tract", "Average_PTSDF", "IBL_Rank"))

prox_heavy_traff <- read_csv(here("spatial_data", "env_health_data", 
                                  "Proximity_to_Heavy_Traffic_Roadways.csv"), 
                             name_repair = janitor::make_clean_names) %>% 
  select(census_tract, proximity_to_heavy_traffic_roadways)
#     col_names=c("County_Name", "Census_Tract", "Prox_Heavy_Traffic_Roadways", "IBL_Rank"))

prox_superfund <- read_csv(here("spatial_data", "env_health_data", "Proximity_to_National_Priorities_List_Facilities_(Superfund_Sites).csv"), 
                           name_repair = janitor::make_clean_names) %>% 
  select( census_tract, average_pnpl)
# col_names=c("County_Name", "Census_Tract", "Average_PNPL", "IBL_Rank"))

toxic_release <- read_csv(here("spatial_data", "env_health_data", "Toxic_Releases_from_Facilities_(RSEI_Model).csv"), 
                          name_repair = janitor::make_clean_names) %>% 
  select(census_tract, average_rsei_concentrations)
# col_names=c("County_Name", "Census_Tract", "Average_RSEI_Concentrations", "IBL_Rank"))

wastewater <- read_csv(here("spatial_data", "env_health_data", "Wastewater_Discharge.csv"), 
                       name_repair = janitor::make_clean_names)  %>% 
  select(census_tract, average_pwdis)
# col_names=c("County_Name", "Census_Tract", "Average_PWDIS", "IBL_Rank"))

diesel_nox <- read_csv(here("spatial_data", "env_health_data", "Diesel_Emission_Levels_of_NOx_(Annual_Tons_Km2).csv"), 
                       name_repair = janitor::make_clean_names) %>% 
  select(census_tract, annual_tons_km2)
# col_names=c("County_Name", "Census_Tract", "Annual Tons Km2 Diesel NOx", "IBL_Rank"))

prox_rmp <- read_csv(here("spatial_data", "env_health_data", "Proximity_to_Risk_Management_Plan_(RMP)_Facililties.csv"), 
                     name_repair = janitor::make_clean_names) %>% 
  select(census_tract, average_prmp)
#  col_names=c("County_Name", "Census_Tract", "Average PRMP", "IBL_Rank"))

ozone <- read_csv(here("spatial_data", "env_health_data", "Ozone_Concentration.csv"), 
                  name_repair = janitor::make_clean_names) %>% 
  select(census_tract, average_ozone_concentration_ppm_km2)
# col_names=c("County_Name", "Census_Tract", "Average Ozone Concentration ppb km2", "IBL_Rank"))


## combine all into one DF 

# first list data frames together 
listdf <- list(lead, pm25, prox_haz_waste, prox_heavy_traff, prox_superfund,
               toxic_release, wastewater, diesel_nox, prox_rmp, ozone)

# remove ranking number 
#listdf <- lapply(listdf, subset, select = -county_name)

# combine data by GEOID 
envdat <- listdf %>% reduce(full_join, by = "census_tract")

# remove first row
envdat <- tail(envdat, -1)

#get only our counties of interest 
envdat <- envdat %>% dplyr::filter(substr(census_tract, 1, 5) 
                                   %in% c("53033", "53053", # king pierce 
                                          "53061", "53035", # snoho kitsap
                                          "53029", "53057")) # island skagit


# rename census tract to geoid and rename pollution columns to be a bit more concise 
envdat <- rename(envdat,c("GEOID" = "census_tract", 
                          "units_with_lead" = "number_units_w_lead_hazard_estimated",
                          "pm2.5_count" = "count",
                          "avg_ptsdf" = "average_ptsdf",
                          "prox_heavy_traffic" = "proximity_to_heavy_traffic_roadways",
                          "avg_pnpl" = "average_pnpl",
                          "avg_rsei" = "average_rsei_concentrations",
                          "avg_pwdis" = "average_pwdis",
                          "diesel_tons_km2" = "annual_tons_km2",
                          "avg_prmp" = "average_prmp",
                          "avg_ozone" = "average_ozone_concentration_ppm_km2")) 

## create polygons for env health data based on geoids 
# we'll use the income data to do this 

# read in income data
tractsKP <- st_read(here("spatial_data", "w_wa_urban_med_income.shp")) %>% 
  st_transform(crs = "EPSG:32610")

st_crs(tractsKP)

# merge env health data to polygons 

# join based on GEOID and put in projection we want for analyses
envdatSP <- merge(tractsKP, envdat, by = "GEOID") %>% 
  st_as_sf() %>%
  st_transform(crs = "EPSG:32610")
  

st_crs(envdatSP)
mapview(envdatSP)

# write shapefile for our env health metrics 
st_write(envdatSP, here("spatial_data", "covariates", "env_health_all_w_wa.shp"), append = FALSE)



##### read in latrine data so we have our points for the buffers 
# same proj 

sites <- st_read(here("clean_data", "latrines.shp")) %>% 
  st_transform(crs = "EPSG:32610")



envdatSP

# plot one of our variables to check 
mapview(list(sites, envdatSP),
        zcol = list(NULL, "avg_pnpl"))

ggplot(data = envdatSP) +
  geom_sf(aes(fill = units_with_lead)) +
  theme_minimal() +
  labs(fill = "number units with lead hazard")

# looks good 

##### now calculate site valeus 

################
# site values for lead 

# rasterize the data 
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
lead_rast <- terra::rasterize(vect(envdatSP), template, field = "units_with_lead")

# initialize new col for data
lead_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(units_lead = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(lead_rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  lead_rast_cropped <- crop(lead_rast, buff_terra)
  lead_rast_masked <- terra::mask(lead_rast_cropped, buff_terra)
  
  # Extract values within the buffer and compute the mean
  extracted_value <- exact_extract(lead_rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  lead_dat$units_lead[i] <- extracted_value
  
  # Convert the raster to a data frame for ggplot2
  lead_rast_df <- as.data.frame(lead_rast_masked, xy = TRUE)
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(lead_rast_cropped, mapping = aes(x = x, y = y, fill = units_with_lead)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(lead_rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  
  print(i)
}

Sys.time()

# Check the final results
print(lead_dat)

# save 
write_csv(lead_dat, here("spatial_data", "covariates", "lead_dat.csv"))







##############  
#site values for pm2.5 


# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "pm2.5_count")

# initialize new col for data
pm2.5_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(pm2.5_count = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  pm2.5_dat$pm2.5_count[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = pm2.5_count)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(pm2.5_dat)

# save
write_csv(pm2.5_dat, here("spatial_data", "covariates", "pm2.5_dat.csv"))





##############  
### site values for PTSDF 

# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "avg_ptsdf")

# initialize new col for data
ptsdf_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(avg_ptsdf = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  ptsdf_dat$avg_ptsdf[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = avg_ptsdf)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(ptsdf_dat)

# write 
write_csv(ptsdf_dat, here("spatial_data", "covariates", "ptsdf_dat.csv"))


##############  
### site values for heavy traffic

# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "prox_heavy_traffic")

# initialize new col for data
traffic_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(prox_heavy_traffic = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  traffic_dat$prox_heavy_traffic[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = prox_heavy_traffic)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(traffic_dat)


# save 
write_csv(traffic_dat, here("spatial_data", "covariates", "traffic_dat.csv"))





##############  
### site values for superfund sites 

# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "avg_pnpl")

# initialize new col for data
pnpl_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(avg_pnpl = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  pnpl_dat$avg_pnpl[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = avg_pnpl)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(pnpl_dat)

# save 
write_csv(pnpl_dat, here("spatial_data", "covariates", "pnpl_dat.csv"))


##############  
### site values for toxic release facilities (rsei)

# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "avg_rsei")

# initialize new col for data
rsei_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(avg_rsei = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  rsei_dat$avg_rsei[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = avg_rsei)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(rsei_dat)

# save 
write_csv(rsei_dat, here("spatial_data", "covariates", "rsei_dat.csv"))


##############  
### site values for wastewater discharge


# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "avg_pwdis")

# initialize new col for data
pwdis_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(avg_pwdis = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  pwdis_dat$avg_pwdis[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = avg_pwdis)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(pwdis_dat)

write_csv(pwdis_dat, here("spatial_data", "covariates", "pwdis_dat.csv"))




##############  
### site values for diesel pollution

# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "diesel_tons_km2")

# initialize new col for data
diesel_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(diesel_tons_km2 = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  diesel_dat$diesel_tons_km2[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = diesel_tons_km2)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(diesel_dat)

#save 
write_csv(diesel_dat, here("spatial_data", "covariates", "diesel_dat.csv"))




##############  
### site values for prmp

# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "avg_prmp")

# initialize new col for data
prmp_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(avg_prmp = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  prmp_dat$avg_prmp[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = avg_prmp)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(prmp_dat)

# save
write_csv(prmp_dat, here("spatial_data", "covariates", "prmp_dat.csv"))





##############  
### site values for ozone

# rasterize the data
template <- rast(ext(envdatSP), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(envdatSP), template, field = "avg_ozone")

# initialize new col for data
ozone_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(avg_ozone = NA_real_)

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, crs(rast)))
  buff_terra <- terra::buffer(buff, width = buffer_radius) 
  
  # need to give it a bit extra extent 
  buff_ext <- ext(buff_terra) 
  buff_ext <- ext(c(
    xmin = buff_ext$xmin - 1000,
    xmax = buff_ext$xmax + 1000,
    ymin = buff_ext$ymin - 1000,
    ymax = buff_ext$ymax + 1000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  ozone_dat$avg_ozone[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = avg_ozone)) +
    #   scale_fill_manual() +
    geom_sf(data = st_as_sf(buff_terra), fill = NA, color = "red", size = 1) +
    geom_sf(data = sites[i,], color = "blue", size = 3) +
    coord_sf(crs = st_crs(rast), datum = st_crs(rast)) +
    labs(title = paste("Site:", i),
         x = "Easting (m)", y = "Northing (m)",
         fill = "Raster Value") +
    theme_minimal()
  
  # Save the plot
  # ggsave(filename = paste0("plots/plot_site_", i, ".png"), plot = p, width = 8, height = 6)
  print(p) 
  Sys.sleep(2)
  print(i)
}

Sys.time()

# Check the final results
print(ozone_dat)


# save
write_csv(ozone_dat, here("spatial_data", "covariates", "ozone_dat.csv"))



#################
# Combine all env health variables and create burden percentiles 

# below code adapted from cesar estien

### bring in data 
lead_dat <- read_csv(here("spatial_data", "covariates", "lead_dat.csv"))
pm2.5_dat <- read_csv(here("spatial_data", "covariates", "pm2.5_dat.csv"))
ptsdf_dat <- read_csv(here("spatial_data", "covariates", "ptsdf_dat.csv"))
traffic_dat <- read_csv(here("spatial_data", "covariates", "traffic_dat.csv"))
pnpl_dat <- read_csv(here("spatial_data", "covariates", "pnpl_dat.csv"))
rsei_dat <- read_csv(here("spatial_data", "covariates", "rsei_dat.csv"))
pwdis_dat <- read_csv(here("spatial_data", "covariates", "pwdis_dat.csv"))
diesel_dat  <- read_csv(here("spatial_data", "covariates", "diesel_dat.csv"))
prmp_dat <- read_csv(here("spatial_data", "covariates", "prmp_dat.csv"))
ozone_dat <- read_csv(here("spatial_data", "covariates", "ozone_dat.csv"))

#now lets make one big pollutant dataframe for the zonal stats
head(lead_dat)
# first list data frames together 
listdfs <- list(lead_dat, pm2.5_dat, ptsdf_dat, traffic_dat, 
                pnpl_dat, rsei_dat, pwdis_dat, 
                diesel_dat, prmp_dat, ozone_dat)

head(lead_dat)

# define a function to keep only desired columns
keep_columns <- function(df) {
  df %>%  # Keep city, site, and remove all other columns
    { if ("scat_id" %in% colnames(.)) dplyr::select(., -scat_id) else . } %>%   
    { if ("prevs_d" %in% colnames(.)) dplyr::select(., -prevs_d) else . } %>%  
    { if ("collctr" %in% colnames(.)) dplyr::select(., -collctr) else . } %>%  
    { if ("month" %in% colnames(.)) dplyr::select(., -month) else . } %>%  
    { if ("date" %in% colnames(.)) dplyr::select(., -date) else . } %>%  
    { if ("certnty" %in% colnames(.)) dplyr::select(., -certnty) else . } %>%  
    { if ("frshnss" %in% colnames(.)) dplyr::select(., -frshnss) else . } %>%  
    { if ("weather" %in% colnames(.)) dplyr::select(., -weather) else . } %>%  
    { if ("notes" %in% colnames(.)) dplyr::select(., -notes) else . } %>%  
    { if ("anl_jll" %in% colnames(.)) dplyr::select(., -anl_jll) else . } %>%  
    { if ("conditn" %in% colnames(.)) dplyr::select(., -conditn) else . } 
}


# apply to all dfs in list 
cleaned_dfs <- purrr::map(listdfs, keep_columns)

# combine data by city/site
pollutants <- cleaned_dfs  %>% reduce(left_join, by = c("city", "site")) 


#center and scale each metric. we'll use this to create the pollution burden score later

# scaling 
pollutants <- pollutants %>%
  mutate(LEAD_PCT = base::scale(units_lead,
                          center=min(units_lead,na.rm=TRUE),
                          scale=diff(range(units_lead,na.rm=TRUE)))[,1]) %>%
  mutate(PM25_PCT = base::scale(pm2.5_count,
                          center=min(pm2.5_count,na.rm=TRUE),
                          scale=diff(range(pm2.5_count,na.rm=TRUE)))[,1]) %>%
  mutate(PTSDF_PCT = base::scale(avg_ptsdf,
                           center=min(avg_ptsdf,na.rm=TRUE),
                           scale=diff(range(avg_ptsdf,na.rm=TRUE)))[,1]) %>%
  mutate(TRAFFIC_PCT = base::scale(prox_heavy_traffic,
                             center=min(prox_heavy_traffic,na.rm=TRUE),
                             scale=diff(range(prox_heavy_traffic,na.rm=TRUE)))[,1]) %>%
  mutate(PNPL_PCT = base::scale(avg_pnpl,
                          center=min(avg_pnpl,na.rm=TRUE),
                          scale=diff(range(avg_pnpl,na.rm=TRUE)))[,1]) %>%
  mutate(TOXIC_PCT = base::scale(avg_rsei,
                           center=min(avg_rsei,na.rm=TRUE),
                           scale=diff(range(avg_rsei,na.rm=TRUE)))[,1]) %>%
  mutate(PWDIS_PCT = base::scale(avg_pwdis,
                           center=min(avg_pwdis,na.rm=TRUE),
                           scale=diff(range(avg_pwdis,na.rm=TRUE)))[,1]) %>%
  mutate(OZONE_PCT = base::scale(avg_ozone,
                           center=min(avg_ozone,na.rm=TRUE),
                           scale=diff(range(avg_ozone,na.rm=TRUE)))[,1]) %>%
  mutate(DIESEL_PCT = base::scale(diesel_tons_km2,
                            center=min(diesel_tons_km2,na.rm=TRUE),
                            scale=diff(range(diesel_tons_km2,na.rm=TRUE)))[,1]) %>%
  mutate(PRMP_PCT = base::scale(avg_prmp,
                          center=min(avg_prmp,na.rm=TRUE),
                          scale=diff(range(avg_prmp,na.rm=TRUE)))[,1]) 

  
#first we need to average the exposures (pm25, diesel, toxic releases (RESI), ozone)
# we'll skip prox traffic because it's too similar to the others

# and environmental risks (superfund sites, risk management plan facilities, hazardous waste facilities, 
# lead risk, wastewater)


pollutants <- pollutants %>%
  mutate(EXPOSURE = (PM25_PCT + DIESEL_PCT +  TOXIC_PCT + OZONE_PCT) / 4) %>%
  mutate(RISK = (LEAD_PCT + PRMP_PCT + PTSDF_PCT + PNPL_PCT + PWDIS_PCT) / 4) 

# write pollutant data 
write_csv(pollutants, here("spatial_data", "covariates", "pollutants_all.csv"))



############
# bring in urbanization covariate data and merge 

urbcovs <- read_csv(here("spatial_data", "covariates", "housing_dat.csv"))
glimpse(pollutants)

pollutants <- read_csv(here("spatial_data", "covariates", "pollutants_all.csv"))

all_covs <- left_join(urbcovs, pollutants, by = c("city", "site")) %>%
  select(-c(1:5,8:13))

head(all_covs)

write_csv(all_covs, here("spatial_data", "covariates", "ALL_ENV_URB_SITES_2000m.csv"))