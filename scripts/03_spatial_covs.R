
# read in spatial data 

# library(pacman)
pacman::p_load(ggplot2, tidyr, dplyr, mapview, lme4, magrittr, lintr, sf, 
               raster, viridis, cowplot, markdown, sf, here, tidycensus,
               crsuggest, terra, spatialEco, readr, ggfortify,  usmap,
               rnaturalearth, rnaturalearthdata, maps, tools, stringr,
               rmapshaper, cowplot, ggrepel, ggspatial, extrafont,exactextractr)

pacman::p_load(dplyr,here,sf,devtools,mapview,here,remotes,tidycensus,crsuggest,
       tidyr,terra,spatialEco,readr,ggfortify,exactextractr,
       spatstat,spdep,landscapemetrics,tmap,viridis,gghighlight,ggcorrplot,
       GGally,magrittr,RStoolbox)

pacman::p_load(here,sf,mapview,remotes,tidycensus,
       tidyr,terra,spatialEco,readr,ggfortify,exactextractr,spatstat,
       spdep,landscapemetrics,tmap,viridis,gghighlight,purrr,plainview,tidyterra,
       crsuggest,janitor,dplyr)

# load in points 
scats <- read_csv(here("data_from_ellie", "field_site_data.csv"),
                      col_names = TRUE) 

colnames(scats) <- (c("scat_id", "previous_id", "collector", "date", "month", "site", "city",
    "lat", "long", "certainty", 
    "freshness", "weather", "condition",
    "notes", "anal_jelly"))

scats <- scats %>% filter(!is.na(scat_id))

# check class of lat long columns
sapply(scats$lat, class)

# turn into sf 
scats <- st_as_sf(scats, coords = c("long", "lat"), crs = 
                    "EPSG:4326")
?st_cast
mapview(scats)

# use utm for geometry
scats <- scats %>%
  st_transform("+proj=utm +zone=10 +datum=WGS84 +units=m")

mapview(scats)


# now let's keep only 1 value for each latrine so we can more easily do the 
# spatial calculations 


sites <- scats %>% distinct(site, .keep_all = TRUE) 

# let's write this to a csv so we can find it later 

st_write(sites, here("clean_data", "latrines."))



########################################
# anthropogenic covariates


###### housing density
# data from: Silvis lab https://silvis.forest.wisc.edu/data/housing-block-change/
# Using 2020 data

wa_housing <- st_read(here("spatial_data",  
                           "WA_block20_change_1990_2020_PLA4.shp"))

# filter to only counties we need
colnames(wa_housing)
wa_housing <- wa_housing %>% dplyr::filter(substr(BLK20, 1, 5) 
                                           %in% c("53033", #king
                                                  "53029", #island
                                                  "53057")) #skagit

# we only need 2020 housing data - select relevant info
wa_housing <- wa_housing %>% dplyr::select(BLK20, WATER20, POP2020,
                                           POPDEN2020, HUDEN2020,
                                           Shape_Leng:geometry)

# crop to actual study area 
sf_use_s2(FALSE)
wa_housing <- st_transform(wa_housing, crs=4326)
wa_housing <- st_crop(wa_housing, c(xmin= -121.7, ymin = 46.7, xmax = -122.8, ymax = 48.5))

#st_write(wa_housing, here("spatial_data", "covariates", "wa_urban_huden_2020.shp"),
  #       append = FALSE)

# let's make all polygons with water (WATER20) NA so they don't get 
# counted in the calculation (otherwise will show up at pop den 0)
wa_housing <- wa_housing %>%
  mutate(POPDEN2020 = ifelse(WATER20 == 1, NA, POPDEN2020))

#mapview(wa_housing)

## now rasterize & calculate site housing density 

# rasterize the data
wa_housing <- st_transform(wa_housing, crs = "EPSG:32610")  
template <- rast(ext(wa_housing), resolution=100, crs="EPSG:32610")
rast <- terra::rasterize(vect(wa_housing), template, field = "POPDEN2020")

plot(rast)
# initialize new col for data
housing_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(POPDEN2020 = NA_real_)

mapview(rast_cropped)
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
    xmin = buff_ext$xmin - 2000,
    xmax = buff_ext$xmax + 2000,
    ymin = buff_ext$ymin - 2000,
    ymax = buff_ext$ymax + 2000))
  buff_terra <- crop(buff_terra, buff_ext)
  
  # crop and mask the raster 
  rast_cropped <- crop(rast, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  housing_dat$POPDEN2020[i] <- extracted_value
  
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, 
                mapping = aes(x = x, y = y, fill = POPDEN2020)) +
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



#write_csv(housing_dat, here("spatial_data", "covariates", "housing_dat.csv"))

########## impervious surface

# read in data
is <-rast(here("spatial_data", "Annual_NLCD_ImpDsc_2022_CU_C1V0.tif"))
st_crs(is)
# first crop to state of washington extent in correct coordinates

# define extent
wa_extent <- ext(c(xmin = -124.848974, xmax = -116.916944, 
               ymin = 45.543541, ymax = 49.002494))  

# create a polygon from the extent 
wa_polygon <- vect(wa_extent, crs = "EPSG:4326")  
?ext

# match raster crs 
wa_extent_albers <- st_transform(st_as_sf(wa_polygon), st_crs(is))

#crop the raster and reproject/crop again to our sf 
is <- crop(is, ext(wa_extent_albers)) %>%
terra::project("EPSG:32610") %>% 
  crop(ext(wa_housing)) %>%       
  mask(ext(wa_housing)) 

plot(is)

# save the cropped raster so we dont have to do this again
#writeRaster(is, here("spatial_data", "washington_NLCD_imp_32610.tif"),
 #           overwrite=TRUE)

is <-rast(here("spatial_data", "washington_NLCD_imp_32610.tif"))


# reproject our SFs
is <- project(is, "EPSG:32610")
wa_housing <- st_transform(wa_housing, crs = st_crs(is)) 
sites <- st_transform(sites, crs = st_crs(is))

is_dat <- st_drop_geometry(sites) %>%
  dplyr::mutate(Impervious = NA_real_)
names(is)
print(is)

summary(values(is))

# set buffer radius
buffer_radius <- 2000

Sys.time()

for (i in 1:nrow(sites)) {
  pt <- sites[i, ]  # iterate through sites 
  
  # create buffer
  buff <- vect(st_transform(pt, st_crs(is)))
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
  rast_cropped <- crop(is, buff_terra)
  rast_masked <- terra::mask(rast_cropped, buff_terra)
  
  # xxtract values within the buffer, calculate mean
  extracted_value <- exact_extract(rast_masked, st_as_sf(buff_terra), fun = "mean", 
                                   weights = "area")
  is_dat$Impervious[i] <- extracted_value
  
  # Plot the results
  
  p <- ggplot() +
    geom_raster(rast_cropped, mapping = aes(x = x, y = y, fill = Annual_NLCD_ImpDsc_2022_CU_C1V0)) +
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
print(is_dat)
print(sites)

write_csv(is_dat, here("data", "covariates", "is_dat.csv"))

