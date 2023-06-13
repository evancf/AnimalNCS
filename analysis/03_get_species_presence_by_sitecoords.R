# Want to get bird and mammal species lists at each site.
# Will do this outside of this code, as it requires IUCN bird and mammal range
# shapefiles that have access controlled by IUCN.

ipak(c("raster", "rgdal", "sf", "tictoc", "beepr"))


# Load mammal range data
mamm_ranges <- shapefile("~/Dropbox/*Science/*Research/*SESYNC/1 Predicting interactions/Data/distribution/range maps/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp")
# This takes about 5 minutes?

# Load bird range data
tic()
bird_ranges <- st_read(dsn = "~/Dropbox/*Science/*Research/*SESYNC/1 Predicting interactions/Data/distribution/range maps/BOTW/BOTW.gdb", layer="All_Species")
toc()
beep(3) # Takes about 18 minutes I think?

# Bring in study coordinates
coord_dat <- read.csv("./data/coord_dat.csv", row.names = 1)

# projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#
# coord_dat <- st_as_sf(coord_dat,
#                         coords = c("long_dec", "lat_dec"),
#                         crs = projcrs)

# Working on projection
mollcrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coord_dat <- coord_dat %>% st_as_sf(coords = c("long_dec", "lat_dec"),
                                    crs = projcrs) %>%
  st_transform(crs = mollcrs)

coord_dat[1,]

# Want a 5km buffer (using appropriate projection for the study site)
coord_buff_500 <- st_buffer(coord_dat, dist = 5000)

# Spatial join and spread to make a big matrix of species presence
mamm_ranges_sf <- st_as_sf(mamm_ranges) %>%
  st_as_sf(crs = projcrs) %>%
  st_transform(crs = mollcrs)

mamm_presence_by_sitecoords <- coord_buff_500 %>% st_join(dplyr::filter(mamm_ranges_sf,
                                                              presence == 1))

mamm_presence_by_sitecoords_wide <- mamm_presence_by_sitecoords[,c("coord_id", "binomial")]
st_geometry(mamm_presence_by_sitecoords_wide) <- NULL
mamm_presence_by_sitecoords_wide <- mamm_presence_by_sitecoords_wide[!duplicated(mamm_presence_by_sitecoords_wide),]
mamm_presence_by_sitecoords_wide$value <- T
mamm_presence_by_sitecoords_wide <- tidyr::spread(data = mamm_presence_by_sitecoords_wide, key = coord_id, value = value, fill = F)

save(mamm_presence_by_sitecoords_wide, file = "./data/spatial data/mamm_presence_by_sitecoords_wide.RData")
Sys.time()


# Spatial join and spread to make a big matrix of species presence
tic()
bird_ranges_sf <- st_as_sf(bird_ranges) %>%
  st_as_sf(crs = projcrs) %>%
  st_transform(crs = mollcrs)

bird_ranges_sf <- dplyr::filter(bird_ranges_sf,
                                PRESENCE == 1) %>% st_zm()

# There are some multisurfaces to deal with
ipak("gdalUtilities")

ensure_multipolygons <- function(X) {
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  st_write(X, tmp1)
  ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- st_read(tmp2)
  st_sf(st_drop_geometry(X), geom = st_geometry(Y))
}

bird_ranges_sf[which(st_geometry_type(bird_ranges_sf) != "MULTIPOLYGON"),
               ] <- bird_ranges_sf[which(st_geometry_type(bird_ranges_sf) != "MULTIPOLYGON"),
                                   ] %>% ensure_multipolygons() %>%
  st_set_crs(st_crs(bird_ranges_sf))


bird_presence_by_sitecoords <- coord_buff_500 %>% st_join(bird_ranges_sf)
bird_presence_by_sitecoords_wide <- bird_presence_by_sitecoords[,c("coord_id", "SCINAME")]
st_geometry(bird_presence_by_sitecoords_wide) <- NULL
bird_presence_by_sitecoords_wide <- bird_presence_by_sitecoords_wide[!duplicated(bird_presence_by_sitecoords_wide),]
bird_presence_by_sitecoords_wide$value <- T
bird_presence_by_sitecoords_wide <- tidyr::spread(data = bird_presence_by_sitecoords_wide, key = coord_id, value = value, fill = F)

save(bird_presence_by_sitecoords_wide, file = "./data/spatial data/bird_presence_by_sitecoords_wide.RData")
toc()

