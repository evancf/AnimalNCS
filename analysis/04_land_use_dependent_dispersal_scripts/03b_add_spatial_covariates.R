# Load packages ----------------------------------------------------------------
sapply(list.files("./R", full.names = T), source)

ipak(c("tidyverse",
       "sf",
       "raster",
       "ggplot2"))

# Load move dataframe and make sf object
move_df <- read.csv("./data/land_use_dependent_dispersal_data/tidy/move_df_with_gee_covs.csv")

mollcrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
move_sf <- move_df %>% st_as_sf(coords = c("location.long", "location.lat"),
                        crs = projcrs) %>%
  st_transform(crs = mollcrs)

move_sf[1,]


# Load ocean and land footprint index
ocean_fp <- raster("./data/spatial data/cumulative_impact_2013.tif") # https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF12B8WBS
#land_ghm <- raster("./data/spatial data/lulc-human-modification-terrestrial-systems_mollweide.tif")
land_fp <- raster("./data/spatial data/HFP2009.tif") # https://datadryad.org/stash/dataset/doi:10.5061/dryad.052q5


# Now extract
land_fp_vals <- raster::extract(land_fp,
                                 move_sf)
ocean_fp_vals <- raster::extract(ocean_fp,
                                 move_sf)

# Bind the footprint columns
land_fp_df <- tibble(land_fp = land_fp_vals)
ocean_fp_df <- tibble(ocean_fp = ocean_fp_vals)

# In order to save this again, we will actually load in the
# original, non-sf version of the movement data
move_df <- move_df %>%
  bind_cols(land_fp_df) %>%
  bind_cols(ocean_fp_df)


# Write this to CSV
write.csv(move_df, file = "./data/land_use_dependent_dispersal_data/tidy/move_df_spat.csv", row.names = F)


