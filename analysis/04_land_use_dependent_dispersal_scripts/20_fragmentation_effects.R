# Packages
devtools::source_gist("71d758f65261a72ab7dc") # gist for ipak

ipak(c("tidyverse",
       "sf"))

# Load Timmers et al. data # https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/fee.2485
# Link to dryad data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.wpzgmsbn9
timm_dat <- read.csv("./data/land_use_dependent_dispersal_data/fragmentation/masterfile_timmers2021_commas_removed.csv",
                     header = T, sep = ";", fill = T)

timm_dat$coord_id <- paste(timm_dat$latitudeN, timm_dat$longitudeE) %>%
  as.factor() %>%
  as.numeric() %>%
  paste0("s", .)

# Pull the coordinates
timm_coords <- timm_dat[,c("coord_id", "latitudeN", "longitudeE")] %>% unique()


# Map these
library(leaflet)
leaflet(timm_coords) %>%
  addProviderTiles('Esri.WorldImagery') %>%
  addMarkers(lng = ~ longitudeE, lat = ~latitudeN)

# Get percent forest cover - because there is a broad range of dates centered
# sort of around 2000, we will just use the year 2000 tree cover from Hansen et al.
mollcrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
timm_sf <- timm_coords %>%
  st_as_sf(coords = c("longitudeE", "latitudeN"),
           crs = projcrs) %>%
  st_transform(crs = mollcrs)

timm_sf[1,]

# Will use a 500 m buffer to roughly match the spatial scale that the
# human footprint data are available for
timm_buff <- st_buffer(timm_sf, dist = 500)

library("rgee")
ee_Initialize(gcs = T, drive = T)

gfc <- ee$Image("UMD/hansen/global_forest_change_2021_v1_9")$select("treecover2000")

# Get the values within the buffer region
gfc_df <- ee_extract(x = gfc,
                          y = timm_buff[, "coord_id"],
                          fun = ee$Reducer$mean(),
                          scale = 30.92, # This is the nominal scale of this dataset
                          sf = T)

# Join these to the timm_dat
# For ease of analysis, will make tree cover relative to 1
timm_gfc <- timm_dat %>%
  left_join(st_drop_geometry(gfc_df)) %>%
  mutate(treecover2000 = treecover2000/100)


# Now let's relate presence probabilities to tree cover in these fragmented
# forest habitats. Use several random effects, as we aim to isolate the effect
# of tree cover, so will correct for the effect of biome, taxonomy, and any
# non-independence associated with particular studies.
library('lme4')
tc_mod <- glmer(Presabs ~ treecover2000 +
                  (1|WWF_Biome) +
                  (1|Binomial_iucn) +
                  (1|Short_ref),
      family = "binomial",
      data = timm_gfc
      )
tc_mod %>% summary()

curve(plogis(fixef(tc_mod)[1] + fixef(tc_mod)[2] * x), col = "blue",
      xlim = c(0, 1),
      ylim = c(0, 1))

plogis(fixef(tc_mod)[1] + fixef(tc_mod)[2] * 0)
plogis(fixef(tc_mod)[1] + fixef(tc_mod)[2] * 1)

# Here's the relationship
# P(presence) = plogis(-0.99 + 1.30 * treecover)
