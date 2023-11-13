source(list.files("./R", full.names = T))

ipak(c("tidyverse", "ggmap", "leaflet",
       "sf", "raster"))

# Pull in data -----------------------------------------------------------------

# Pull data from GROA.
# Further info and metadata here https://github.com/forc-db/GROA
g_sites <- read.csv("https://raw.githubusercontent.com/forc-db/GROA/master/data/sites.csv")[,-1] %>% tibble()
g_biomass <- read.csv("https://raw.githubusercontent.com/forc-db/GROA/master/data/biomass_litter_CWD.csv") %>% tibble()


# Pull GPFC dataset
# These data are available here: https://zenodo.org/records/6555216
pl_sites <- readxl::read_excel("./data/GPFC_database.xlsx", sheet = 2)
pl_biomass <- readxl::read_excel("./data/GPFC_database.xlsx", sheet = 3)



# Reconcile names --------------------------------------------------------------
g_sites %>% colnames()
g_biomass %>% colnames()

pl_sites %>% colnames()
pl_biomass %>% colnames()

colnames(g_biomass)[colnames(g_biomass) %in% colnames(pl_biomass)]
colnames(pl_biomass)[colnames(pl_biomass) %in% colnames(g_biomass)]

# Let's convert points to underscore in groa
colnames(g_sites) <- gsub(".", "_", colnames(g_sites), fixed = T)
colnames(g_biomass) <- gsub(".", "_", colnames(g_biomass), fixed = T)

# Convert some gpfc colnames to the groa colnames

# Location info
pl_sites <- pl_sites %>%
  rename(long_dec = longitude,
         lat_dec = latitude,
         site_country = country,
         site_state = state)

# Reconcile column names
pl_biomass <- pl_biomass %>%
  rename(mean_ha = agc,
         stand_age = age)



# Convert to a character for the study and plot names, making sure it's clear
# which dataset these came from
pl_biomass <- pl_biomass %>%
  mutate(plot_id = paste0("pl", plot_id),
         site_id = paste0("pl", site_id),
         study_id = paste0("pl", study_id)) %>%
  mutate(monoculture = "yes") # Also want to add a column saying whether it's a monoculture
pl_sites <- pl_sites %>%
  mutate(site_id = paste0("pl", site_id),
         study_id = paste0("pl", study_id))

g_biomass <- g_biomass %>%
  mutate(plot_id = paste0("g", plot_id),
         site_id = paste0("g", site_id),
         study_id = paste0("g", study_id)) %>%
  mutate(monoculture = "no") # Also want to add a column saying whether it's a monoculture
g_sites <- g_sites %>%
  mutate(site_id = paste0("g", site_id),
         study_id = paste0("g", study_id))



# Merge the datasets -----------------------------------------------------------

# Join each individually, then bind rows
pl_dat <- pl_biomass %>% left_join(pl_sites)
g_dat <- g_biomass %>% left_join(g_sites)

all_dat <- bind_rows(g_dat, pl_dat)




# Get site-level values --------------------------------------------------------

# First get a unique value for the coordiantes. This is (very rarely) different
# from the site id
all_dat <- all_dat %>%
  mutate(coord_id = paste(lat_dec, long_dec, sep = "_"))

# Next get a tibble of just the unique coordinates
coord_dat <- all_dat %>%
  dplyr::select(coord_id, long_dec, lat_dec) %>%
  unique()


# This commented out part has some potential code for checking how many
# sites there are versus how many plots, versus how many unique coordinates
# %>%
# mutate(coord_id = paste0(lat_dec, long_dec)) %>%
# dplyr::select(-long_dec, lat_dec)
# table(list(coord_dat$plot_id, coord_dat$coord_id)) %>% table()
# length(unique(coord_dat$plot_id))
# length(unique(coord_dat$site_id))
# length(unique(coord_dat$coord_id))
# table(list(coord_dat$coord_id, coord_dat$site_id)) %>% table()

# Note that there are missing coordinates for various plots/sites. This
# means there is an NA in there in coord_dat. Will filter this away
# so that we can move ahead with analyses that use spatial data
coord_dat <- coord_dat %>%
  filter(complete.cases(coord_dat))


# Will write this to csv for use later
write.csv(file = "./data/coord_dat.csv", coord_dat)


# Working on projection
mollcrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coord_dat <- coord_dat %>% st_as_sf(coords = c("long_dec", "lat_dec"),
                                        crs = projcrs) %>%
  st_transform(crs = mollcrs)

coord_dat[1,]


# Get the ecoregions and realms
# See details in 01_spatial_data_from_gee.R
ecoreg <- read_sf("./data/spatial data/vector/Ecoregions2017/Ecoregions2017.shp")
sf_use_s2(FALSE)
ecoreg_dat <- coord_dat %>%
  st_transform(crs = projcrs) %>%
  st_join(ecoreg)

write.csv(file = "./data/ecoreg_dat.csv",
          ecoreg_dat %>%
            st_drop_geometry() %>%
            arrange(coord_id))


# Get a buffer around the plots
coord_buff_500 <- st_buffer(coord_dat, dist = 500)
coord_buff_2500 <- st_buffer(coord_dat, dist = 2500)
coord_buff_10000 <- st_buffer(coord_dat, dist = 10000)



### Forest cover
library("rgee")
ee_Initialize(gcs = T, drive = T)

gfc <- ee$Image("UMD/hansen/global_forest_change_2021_v1_9")$select("treecover2000")

# Get the values within the buffer region
gfc_df <- ee_extract(x = gfc,
                     y = coord_buff_500[1:1000, "coord_id"],
                     fun = ee$Reducer$mean(),
                     scale = 30.92, # This is the nominal scale of this dataset
                     sf = T)

# Have to do this in chunks for some reason
gfc_df <- gfc_df %>%
  bind_rows(ee_extract(x = gfc,
                       y = coord_buff_500[1001:2000, "coord_id"],
                       fun = ee$Reducer$mean(),
                       scale = 30.92, # This is the nominal scale of this dataset
                       sf = T))
gfc_df <- gfc_df %>%
  bind_rows(ee_extract(x = gfc,
                       y = coord_buff_500[2001:nrow(coord_buff_500), "coord_id"],
                       fun = ee$Reducer$mean(),
                       scale = 30.92, # This is the nominal scale of this dataset
                       sf = T))




### Potential NPP

npp <- ee$ImageCollection("MODIS/006/MOD17A3HGF")$
  select("Npp")$
  filterDate("2001-01-01", "2011-01-01")$
  mean()

# Get the values within the buffer region
npp_10000_max_df <- ee_extract(x = npp,
                               y = coord_buff_10000[1:1000, "coord_id"],
                               fun = ee$Reducer$max(),
                               scale = 500, # This is the nominal scale of this dataset
                               sf = T)

# Have to do this in chunks for some reason
npp_10000_max_df <- npp_10000_max_df %>%
  bind_rows(ee_extract(x = npp,
                       y = coord_buff_10000[1001:2000, "coord_id"],
                       fun = ee$Reducer$max(),
                       scale = 500, # This is the nominal scale of this dataset
                       sf = T))
npp_10000_max_df <- npp_10000_max_df %>%
  bind_rows(ee_extract(x = npp,
                       y = coord_buff_10000[2001:nrow(coord_buff_500), "coord_id"],
                       fun = ee$Reducer$max(),
                       scale = 500, # This is the nominal scale of this dataset
                       sf = T))

npp_10000_max_df <- npp_10000_max_df %>%
  rename("npp_10000_max_average" = "Npp")



### Human modification

ghm <- ee$ImageCollection("CSP/HM/GlobalHumanModification")

# Get the values within the buffer region
ghm_2500_mean_df <- ee_extract(x = ghm,
                               y = coord_buff_2500[1:1000, "coord_id"],
                               fun = ee$Reducer$mean(),
                               scale = 1000, # This is the nominal scale of this dataset
                               sf = T)

# Have to do this in chunks for some reason
ghm_2500_mean_df <- ghm_2500_mean_df %>%
  bind_rows(ee_extract(x = ghm,
                       y = coord_buff_2500[1001:2000, "coord_id"],
                       fun = ee$Reducer$mean(),
                       scale = 1000, # This is the nominal scale of this dataset
                       sf = T))
ghm_2500_mean_df <- ghm_2500_mean_df %>%
  bind_rows(ee_extract(x = ghm,
                       y = coord_buff_2500[2001:nrow(coord_buff_500), "coord_id"],
                       fun = ee$Reducer$mean(),
                       scale = 1000, # This is the nominal scale of this dataset
                       sf = T))



### Land footprint index
land_fp <- raster("./data/spatial data/raster/HFP2009.tif") # https://datadryad.org/stash/dataset/doi:10.5061/dryad.052q5

land_fp_vals <- raster::extract(land_fp,
                                coord_dat)


### Temp and precip

### Temperature

temp <- ee$Image("WORLDCLIM/V1/BIO")$select("bio01")

# Get the values within the buffer region
temp_500_mean_df <- ee_extract(x = temp,
                               y = coord_buff_500[1:1000, "coord_id"],
                               fun = ee$Reducer$mean(),
                               scale = 927.67, # This is the nominal scale of this dataset
                               sf = T)

# Have to do this in chunks for some reason
temp_500_mean_df <- temp_500_mean_df %>%
  bind_rows(ee_extract(x = temp,
                       y = coord_buff_500[1001:2000, "coord_id"],
                       fun = ee$Reducer$mean(),
                       scale = 927.67, # This is the nominal scale of this dataset
                       sf = T))
temp_500_mean_df <- temp_500_mean_df %>%
  bind_rows(ee_extract(x = temp,
                       y = coord_buff_500[2001:nrow(coord_buff_500), "coord_id"],
                       fun = ee$Reducer$mean(),
                       scale = 927.67, # This is the nominal scale of this dataset
                       sf = T))

temp_500_mean_df <- temp_500_mean_df %>%
  dplyr::rename(mat_C = bio01) %>%
  mutate(mat_C = mat_C / 10) # For some reason this is in units of 0.1 degree C

### Precipitation

precip <- ee$Image("WORLDCLIM/V1/BIO")$select("bio12")

# Get the values within the buffer region
precip_500_mean_df <- ee_extract(x = precip,
                                 y = coord_buff_500[1:1000, "coord_id"],
                                 fun = ee$Reducer$mean(),
                                 scale = 927.67, # This is the nominal scale of this dataset
                                 sf = T)

# Have to do this in chunks for some reason
precip_500_mean_df <- precip_500_mean_df %>%
  bind_rows(ee_extract(x = precip,
                       y = coord_buff_500[1001:2000, "coord_id"],
                       fun = ee$Reducer$mean(),
                       scale = 927.67, # This is the nominal scale of this dataset
                       sf = T))
precip_500_mean_df <- precip_500_mean_df %>%
  bind_rows(ee_extract(x = precip,
                       y = coord_buff_500[2001:nrow(coord_buff_500), "coord_id"],
                       fun = ee$Reducer$mean(),
                       scale = 927.67, # This is the nominal scale of this dataset
                       sf = T))

precip_500_mean_df <- precip_500_mean_df %>%
  dplyr::rename(tap_mm = bio12)





# Add all the site-level values to dataset -------------------------------------

# First will join all the site-level values to coord_dat
coord_covars <- coord_dat %>%
  st_drop_geometry() %>%
  left_join(st_drop_geometry(gfc_df)) %>% # Forest cover 2000
  #mutate(gpp = gpp_df$gpp) %>% # Potential GPP
  left_join(st_drop_geometry(npp_10000_max_df %>% dplyr::select(-starts_with("X")))) %>% # Maximum NPP
  left_join(st_drop_geometry(ghm_2500_mean_df)) %>% # Human modification
  mutate(fp = land_fp_vals) %>%  # Human footprint values
  left_join(st_drop_geometry(temp_500_mean_df)) %>% # Worldclim temp
  left_join(st_drop_geometry(precip_500_mean_df)) # Worldclim precip

coord_covars <- coord_covars %>%
  dplyr::arrange(coord_id)

write.csv(file = "./data/coord_covars.csv", coord_covars)

# Write out all_dat as well ----------------------------------------------------

write.csv(file = "./data/all_dat.csv", all_dat)
