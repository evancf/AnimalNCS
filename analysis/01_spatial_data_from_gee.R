# Our goal here is to get several global rasters in the same
# resolution / crs saved locally from the GEE catalogue

# Toward the bottom are some other rasters that can be
# downloaded directly from the linked URLs

# Load a couple functions
source(list.files("./R", full.names = T))

topwd <- getwd()

ipak(c("tidyverse","rgee"))
# See README notes about Google Earth Engine
ee_Initialize(gcs = T, drive = T)


# Functions
merge_project_gee_rasters <- function(x){
  t0 <- Sys.time()
  filename <- deparse(substitute(x)) %>%
    gsub("_raster", ".tif", ., fixed = T)

  y <- terra::merge(terra::rast(x[1]),
                    terra::rast(x[2])) %>%
    terra::project(dd_terra)
  terra::writeRaster(y,
                     filename,
                     overwrite = T)
  t1 <- Sys.time()
  print("Done!")
  print(t1-t0)
}


# Will use this as the template raster for the other raster data that is downloaded
setwd(topwd)
dd_terra <- terra::rast("./data/spatial data/raster_outputs/dispersal_disruption2020.tif")


# Want to work in this via_gee folder
setwd(topwd)
setwd("./data/spatial data/raster/via_gee/")

# Global forest watch treecover 2000 -------------------------------------------
if(!"gfc.tif" %in% list.files()){
  gfc_raster <- ee$
    Image("UMD/hansen/global_forest_change_2021_v1_9")$select("treecover2000")$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "gfc.tif")
  # Now save a projected file
  merge_project_gee_rasters(gfc_raster)
}


# Glad loss raster -------------------------------------------------------------
if(!"glad_loss2020.tif" %in% list.files()){
  glad_loss2020_raster <- ee$
    Image("projects/glad/GLCLU2020/Forest_loss")$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "glad_loss2020.tif")
  # Now save a projected file
  merge_project_gee_rasters(glad_loss2020_raster)
}

# Glad gain raster -------------------------------------------------------------
if(!"glad_gain2020.tif" %in% list.files()){
  glad_gain2020_raster <- ee$
    Image("projects/glad/GLCLU2020/Forest_gain")$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "glad_gain2020.tif")
  # Now save a projected file
  merge_project_gee_rasters(glad_gain2020_raster)
}

# Also make a change raster
setwd(topwd)
setwd("./data/spatial data/raster/via_gee/")
if(!"gfc2020.tif" %in% list.files()){
  gfc_stars <- read_stars("gfc.tif", proxy = F)
  glad_loss2020 <- read_stars("glad_loss2020.tif", proxy = F)
  glad_gain2020 <- read_stars("glad_gain2020.tif", proxy = F)

  gfc2020_stars <- gfc_stars/100 - glad_loss2020 + glad_gain2020

  # Not sure if we'll have the memory for this
  gfc2020_stars[[1]][gfc2020_stars[[1]] > 1] <- 1
  gfc2020_stars[[1]][gfc2020_stars[[1]] < 0] <- 0

  write_stars(gfc2020_stars, dsn = "gfc2020.tif")
}


# Temperature raster -------------------------------------------------------------
if(!"amt.tif" %in% list.files()){
  amt_raster <- ee$
    Image("WORLDCLIM/V1/BIO")$select("bio01")$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "amt.tif")
  # Now save a projected file
  merge_project_gee_rasters(amt_raster)
}


# Precipitation raster -------------------------------------------------------------
if(!"tap.tif" %in% list.files()){
  tap_raster <- ee$
    Image("WORLDCLIM/V1/BIO")$select("bio12")$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "tap.tif")
  # Now save a projected file
  merge_project_gee_rasters(tap_raster)
}



# Potential NPP -------------------------------------------
setwd(topwd)
setwd("./data/spatial data/raster/via_gee/")
if(!"npp_max.tif" %in% list.files()){
  npp_max_raster <- ee$
    ImageCollection("MODIS/006/MOD17A3HGF")$
    filterDate("2001-01-01", "2011-01-01")$
    select("Npp")$
    mean()$
    reduceNeighborhood(reducer = ee$Reducer$max(),
                       kernel = ee$Kernel$circle(radius = 10000, units = "meters")) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "npp_max.tif")
  # Now save a projected file
  merge_project_gee_rasters(npp_max_raster)
}




# Will also get some other rasters set up --------------------------------------

# Potential restoration sites from Griscom et al. 2017 PNAS
# https://zenodo.org/record/883444
setwd(topwd)
setwd("./data/spatial data/raster/")
if(!"ncs_restore_sites.tif" %in% list.files()){
  ncs_restore_sites <- terra::rast("./NCS_Refor11_map/NCS_Refor11_map.tif") %>%
    terra::project(dd_terra)
  terra::writeRaster(ncs_restore_sites,
                     "ncs_restore_sites.tif",
                     overwrite = T)
  rm(ncs_restore_sites)
}

# Human footprint
# https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint
setwd(topwd)
setwd("./data/spatial data/raster/")
if(!"fp.tif" %in% list.files()){
  fp <- terra::rast("./HFP2009.tif") %>%
    terra::project(dd_terra)
  terra::writeRaster(fp,
                     "fp.tif",
                     overwrite = T)
  rm(fp)
}



# Some other masks as rasters --------------------------------------------------

# Because we'll need to use the raster package, we will
# get a different template
setwd(topwd)
dd_raster <- raster::raster("./data/spatial data/raster_outputs/dispersal_disruption2020.tif")

# Ecoregions dataset available via this link: https://storage.googleapis.com/teow2016/Ecoregions2017.zip
# Also see the "About" tab here: https://ecoregions.appspot.com/

setwd("./data/spatial data/vector/")
library("sf")
ecoreg <- read_sf("Ecoregions2017/Ecoregions2017.shp")

ecoreg$BIOME_NAME %>% table()

# Create some columns
ecoreg <- ecoreg %>%
  mutate(forest = ifelse(grepl("Forest", BIOME_NAME), 1, 0)) %>%
  mutate(savanna = ifelse(grepl("Savanna", BIOME_NAME), 1, 0)) %>%
  mutate(tropical = ifelse(grepl("Tropical", BIOME_NAME), 1, 0)) %>%
  mutate(for_and_sav = savanna + forest) %>%
  mutate(trop_for = forest * tropical) %>%
  mutate(trop_for_and_sav = for_and_sav * tropical) %>%
  mutate(for_and_sav = ifelse(for_and_sav == 1, 1, NA)) %>%
  mutate(trop_for = ifelse(trop_for == 1, 1, NA)) %>%
  mutate(trop_for_and_sav = ifelse(trop_for_and_sav == 1, 1, NA))

# Transform
proj_crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
ecoreg <- ecoreg %>%
  st_transform(proj_crs)


# Save this as a tif so that it doesn't have to be done again
setwd(topwd)
setwd("./data/spatial data/raster/")
if(!"for_and_sav_biomes.tif" %in% list.files()){
  library("fasterize")
  for_and_sav_biomes <- fasterize(sf = ecoreg,
                                  raster = dd_raster,
                                  field = "for_and_sav")
  raster::writeRaster(for_and_sav_biomes,
                      filename = "for_and_sav_biomes.tif",
                      overwrite = T)
  # Will need to project this to match the proj crs
  for_and_sav_biomes <- terra::rast("for_and_sav_biomes.tif") %>%
    terra::project(dd_terra, method = "near")
  terra::writeRaster(for_and_sav_biomes,
                     "for_and_sav_biomes.tif",
                     overwrite = T)
}


setwd(topwd)
setwd("./data/spatial data/raster/")
if(!"trop_for_and_sav_biomes.tif" %in% list.files()){
  library("fasterize")
  trop_for_biomes <- fasterize(sf = ecoreg,
                               raster = dd_raster,
                               field = "trop_for")
  raster::writeRaster(trop_for_biomes,
                      filename = "trop_for_biomes.tif",
                      overwrite = T)
  # Will need to project this to match the proj crs
  trop_for_biomes <- terra::rast("trop_for_biomes.tif") %>%
    terra::project(dd_terra, method = "near")
  terra::writeRaster(trop_for_biomes,
                     "trop_for_biomes.tif",
                     overwrite = T)
}

setwd(topwd)
setwd("./data/spatial data/raster/")
if(!"trop_for_and_sav_biomes.tif" %in% list.files()){
  library("fasterize")
  trop_for_and_sav_biomes <- fasterize(sf = ecoreg,
                                       raster = dd_raster,
                                       field = "trop_for_and_sav")
  raster::writeRaster(trop_for_and_sav_biomes,
                      filename = "trop_for_and_sav_biomes.tif",
                      overwrite = T)
  # Will need to project this to match the proj crs
  trop_for_and_sav_biomes <- terra::rast("trop_for_and_sav_biomes.tif") %>%
    terra::project(dd_terra, method = "near")
  terra::writeRaster(trop_for_and_sav_biomes,
                     "trop_for_and_sav_biomes.tif",
                     overwrite = T)
}

