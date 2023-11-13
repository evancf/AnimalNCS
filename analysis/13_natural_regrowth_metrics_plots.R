source(list.files("./R", full.names = T))

ipak(c("tidyverse", "stars"))

# Set up working directories
top_wd <- getwd()

data_wd <- paste0(top_wd,"/data/")
outputs_wd <- paste0(top_wd,"/outputs/")
raster_wd <- paste0(top_wd, "/data/spatial data/raster/")
raster_gee_wd <- paste0(top_wd, "/data/spatial data/raster/via_gee/")
raster_outputs_wd <- paste0(top_wd, "/data/spatial data/raster_outputs")
raster_intermediary_wd <- paste0(top_wd, "/data/spatial data/raster_intermediary")

# data_wd <- "/nobackup1/efricke/VertNBS_data"
# outputs_wd <- "/nobackup1/efricke/VertNBS_data"
# raster_wd <- "/nobackup1/efricke/VertNBS_data"
# raster_gee_wd <- "/nobackup1/efricke/VertNBS_data"
# raster_outputs_wd <- "/nobackup1/efricke/VertNBS_visuals"
# raster_intermediary_wd <- "/nobackup1/efricke/VertNBS_data/intermediary_rasters"


# A few objects to help with maps and metrics ----------------------------------
# # # What is the cell size for this projection
# setwd(raster_outputs_wd)
# proj_crs_cell_size <- st_area(read_stars("c30_2020_current_trop_for_and_sav_biomes.tif"), proxy = F)
# range(proj_crs_cell_size$area, na.rm = T)
# # Units: [m^2]
# # [1] 708264.2 708264.2
proj_crs_cell_size <- 708264.2 / 10000

# First get a continents polygon
# Want to match this stars raster
setwd(raster_outputs_wd)
c30_2020_current_trop_for_and_sav_biomes <- read_stars("c30_2020_current_trop_for_and_sav_biomes.tif")
library("rnaturalearth")
sf_use_s2(FALSE)
world <- ne_countries() %>%
  st_as_sf() %>%
  filter(continent != "Antarctica") %>%
  st_union() %>%
  st_transform(crs = st_crs(c30_2020_current_trop_for_and_sav_biomes))
focal_bbox <- st_bbox(c30_2020_current_trop_for_and_sav_biomes)
sf_use_s2(T)

# A version of theme_void that avoids some other formatting issues
theme_blank <- function(){
  theme(line = element_blank(),
        rect = element_blank(),

        axis.line.x = element_blank(),
        axis.line.y = element_blank(),

        axis.title.x = element_blank(),
        axis.title.y = element_blank(),

        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),

        axis.text.x = element_blank(),
        axis.text.y = element_blank())
}

# Set up for hatch marks in savanna biomes
setwd(top_wd)
setwd("./data/spatial data/vector/")
ecoreg <- read_sf("Ecoregions2017/Ecoregions2017.shp")

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
  mutate(trop_for_and_sav = ifelse(trop_for_and_sav == 1, 1, NA)) %>%
  mutate(trop_sav = ifelse(savanna * tropical == 1, 1, NA))

sf_use_s2(FALSE)
trop_sav_sf <- ecoreg %>%
  filter(trop_sav == 1) %>%
  st_union() %>%
  st_simplify(dTolerance = 0.01)
sf_use_s2(TRUE)

trop_sav_poly <- st_cast(trop_sav_sf, "POLYGON")
hatch_list <- list()
library("HatchedPolygons")
for(i in 1:length(trop_sav_poly)){
  # A little work around to avoid errors for very small polygons
  hatch_list[[i]] <- tryCatch(hatched.SpatialPolygons(trop_sav_poly[[i]], density = 0.5), error=function(e){})
}
hatch <- do.call(rbind, hatch_list)
hatch <- hatch %>% sf::st_set_crs(st_crs(trop_sav_sf))
# Transform to match maps
proj_crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
hatch <- hatch %>%
  st_transform(proj_crs)
rm(ecoreg, trop_sav_poly, trop_sav_sf, hatch_list)



# Maps -------------------------------------------------------------------------

# Current regrowth potential
setwd(raster_outputs_wd)
c30_2020_current_trop_for_and_sav_biomes <- read_stars("c30_2020_current_trop_for_and_sav_biomes.tif", proxy = F) + 0

setwd(top_wd)
# pdf(file = "./outputs/figures/regrowth potential 2020.pdf",
#     width = 4.76 * 7.24 / 4.76,
#     height = 1.8 * 7.24 / 4.76)
png(file = "./outputs/figures/regrowth potential 2020.png",
    width = 4.76 * 7.24 / 4.76,
    height = 1.8 * 7.24 / 4.76,
    units = "in", res = 440)

library(colorspace)
ggplot() +
  geom_stars(data = c30_2020_current_trop_for_and_sav_biomes,
             downsample = 20,
             aes(x = x, y = y, fill = c30_2020_current_trop_for_and_sav_biomes.tif)) +
  scale_fill_gradientn(values = c(0, 0.4, 1),
                         colours = c("lightyellow", "forestgreen", "darkgreen"),
                         na.value = "white") +
  guides(fill = guide_colorbar(title="Current regrowth\npotential (Mg/ha/yr)",
                               barwidth = 12,
                               title.vjust = 1.7,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,-0.15),
        legend.direction = "horizontal") +
  theme(plot.margin=grid::unit(c(0,0,5,0), "mm")) +
  ylim(st_bbox(world)[2] * .7, st_bbox(world)[4] * .58) +
  xlim(st_bbox(world)[1] * .9, st_bbox(world)[3] * 1) +
  coord_equal() +
  theme_blank() +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA) +
  geom_sf(data = hatch, color = "white", linewidth = 0.1)
dev.off()

rm(c30_2020_current_trop_for_and_sav_biomes)





# Gains and losses since 2000
setwd(raster_outputs_wd)
list.files()
c30_2000_current_trop_for_and_sav_biomes <- read_stars("c30_2000_current_trop_for_and_sav_biomes.tif", proxy = F)
c30_2020_current_trop_for_and_sav_biomes <- read_stars("c30_2020_current_trop_for_and_sav_biomes.tif", proxy = F)

# Quick diversion to get a sense for the spatial breadth of gains vs losses
diff_ratio <- (c30_2020_current_trop_for_and_sav_biomes[[1]] - c30_2000_current_trop_for_and_sav_biomes[[1]]) / c30_2000_current_trop_for_and_sav_biomes[[1]]
diff_ratio_density <- hist(diff_ratio, breaks = c(-1,-0.05, 0, 0.05, 3))
diff_ratio_density$density[1] / (diff_ratio_density$density[1] + diff_ratio_density$density[4])
rm(diff_ratio)
gc

c30_diff_current_trop_for_and_sav_biomes <- c30_2020_current_trop_for_and_sav_biomes - c30_2000_current_trop_for_and_sav_biomes
rm(c30_2000_current_trop_for_and_sav_biomes)
rm(c30_2020_current_trop_for_and_sav_biomes)


# Also want to calculate the losses in these areas
# Load in the ncs sites
# Consider tropical forest and savanna sites
setwd(raster_wd)
ncs_restore_sites <- read_stars("ncs_restore_sites.tif", proxy = F)

# Will mask to tropical forest and savanna
ncs_restore_sites[[1]] <- ncs_restore_sites[[1]] * read_stars("trop_for_and_sav_biomes.tif", proxy = F)[[1]]

# Want to convert to NA all the pixels that are not potential restoration sites
ncs_restore_sites[[1]][ncs_restore_sites[[1]] == 0] <- NA

#
c30_2020_ncs <- c30_2020_current_trop_for_and_sav_biomes
c30_2020_ncs[[1]] <- c30_2020_ncs[[1]] * ncs_restore_sites[[1]]
c30_2000_ncs <- c30_2000_current_trop_for_and_sav_biomes
c30_2000_ncs[[1]] <- c30_2000_ncs[[1]] * ncs_restore_sites[[1]]
diff_ratio <- (c30_2020_ncs[[1]] - c30_2000_ncs[[1]]) / c30_2000_ncs[[1]]
# Let's see how many pixels with large magnitude change (>+/- 5%)
diff_ratio_density <- hist(diff_ratio, breaks = c(-1,-0.05, 0, 0.05, 1.5))
gc()
mean(diff_ratio, na.rm = T)

c30_diff_ncs <- c30_diff_current_trop_for_and_sav_biomes
c30_diff_ncs[[1]] <- c30_diff_ncs[[1]] * ncs_restore_sites[[1]]
mean(c30_diff_ncs[[1]], na.rm = T)
rm(c30_diff_ncs)


setwd(top_wd)
# pdf(file = "./outputs/figures/gained and lost regrowth potential.pdf",
#     width = 4.76 * 7.24 / 4.76,
#     height = 1.8 * 7.24 / 4.76)
png(file = "./outputs/figures/gained and lost regrowth potential.png",
    width = 4.76 * 7.24 / 4.76,
    height = 1.8 * 7.24 / 4.76,
    units = "in", res = 440)

library(colorspace)
ggplot() +
  geom_stars(data = c30_diff_current_trop_for_and_sav_biomes,
             downsample = 20,
             aes(x = x, y = y, fill = c30_2020_current_trop_for_and_sav_biomes.tif)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,
                                   na.value = "white",
                                   p1 = 0.6,
                                   p2 = 0.6,
                                   p3 = 0.6,
                                   p4 = 0.6) +
  guides(fill = guide_colorbar(title="Change in regrowth potential\nsince 2000 (Mg/ha/yr)",
                               barwidth = 12,
                               title.vjust = 1.7,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(#legend.position="bottom",
    legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,-0.15),
        legend.direction = "horizontal") +
  theme(plot.margin=grid::unit(c(0,0,5,0), "mm")) +
  ylim(st_bbox(world)[2] * .7, st_bbox(world)[4] * .58) +
  xlim(st_bbox(world)[1] * .9, st_bbox(world)[3] * 1) +
  coord_equal() +
  theme_blank() +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA)
dev.off()

rm(c30_diff_current_trop_for_and_sav_biomes)




# Rewilding potential
setwd(raster_outputs_wd)
c30_2020_diff_trop_for_and_sav_biomes <- read_stars("c30_2020_diff_trop_for_and_sav_biomes.tif", proxy = F) + 0
#c30_2020_diff_trop_for_biomes <- read_stars("c30_2020_diff_trop_for_biomes.tif", proxy = F) + 0


setwd(top_wd)
# pdf(file = "./outputs/figures/rewilding potential 2020.pdf",
#     width = 4.76 * 7.24 / 4.76,
#     height = 1.8 * 7.24 / 4.76)
png(file = "./outputs/figures/rewilding potential 2020.png",
    width = 4.76 * 7.24 / 4.76,
    height = 1.8 * 7.24 / 4.76,
    units = "in", res = 440)

library(colorspace)
ggplot() +
  geom_stars(data = c30_2020_diff_trop_for_and_sav_biomes,
             downsample = 20,
             aes(x = x, y = y, fill = c30_2020_diff_trop_for_and_sav_biomes.tif)) +
  # geom_stars(data = c30_2020_diff_trop_for_biomes,
  #            downsample = 20,
  #            aes(x = x, y = y, fill = c30_2020_diff_trop_for_biomes.tif)) +
  scale_fill_viridis_c(option = "plasma",
                       na.value = "white") +
  guides(fill = guide_colorbar(title="Carbon accumulation potential of biodiversity\n& connectivity restoration (Mg/ha/yr)",
                               barwidth = 12,
                               title.vjust = 1.7,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,-0.15),
        legend.direction = "horizontal") +
  theme(plot.margin=grid::unit(c(0,0,5,0), "mm")) +
  ylim(st_bbox(world)[2] * .7, st_bbox(world)[4] * .58) +
  xlim(st_bbox(world)[1] * .9, st_bbox(world)[3] * 1) +
  coord_equal() +
  theme_blank() +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA) +
  geom_sf(data = hatch, color = "white", linewidth = 0.1)
dev.off()

#rm(c30_2020_diff_trop_for_and_sav_biomes)

# Insets?

ggplot() +
  geom_stars(data = c30_2020_diff_trop_for_and_sav_biomes,
             downsample = 2,
             aes(x = x, y = y, fill = c30_2020_diff_trop_for_and_sav_biomes.tif)) +
  # geom_stars(data = c30_2020_diff_trop_for_biomes,
  #            downsample = 2,
  #            aes(x = x, y = y, fill = c30_2020_diff_trop_for_biomes.tif)) +
  # scale_fill_gradientn(values = c(0, 0.3, 1),
  #                      colours = c("white", "lightblue", "darkblue"),
  #                      na.value = "white") +
  scale_fill_viridis_c(option = "magma",
                       na.value = "white") +
  guides(fill = guide_colorbar(title="Regrowth potential\n(Mg/ha/yr)",
                               barwidth = 12,
                               title.vjust = 1.7,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,-0.15),
        legend.direction = "horizontal") +
  theme(plot.margin=grid::unit(c(0,0,5,0), "mm")) +
  ylim(st_bbox(world)[2] * .275, st_bbox(world)[2] * 0) +
  xlim(st_bbox(world)[1] * .9, st_bbox(world)[3] * 1) +
  coord_equal() +
  theme_blank()




# Mapping uncertainty in terms of the magnitude of 95% credible interval ----------
# Uncertainty: current
setwd(raster_outputs_wd)
current_lo <- read_stars("c30_2020_current_trop_for_and_sav_biomes_lo_fullmap.tif",
                         proxy = F) + 0

current_hi <- read_stars("c30_2020_current_trop_for_and_sav_biomes_hi_fullmap.tif",
                         proxy = F) + 0

current_95 <- current_hi
current_95[[1]] <- current_hi[[1]] - current_lo[[1]]

rm(current_lo, current_hi)
gc()

# Want to plot pixels without vales  as white not grey
current_95[[1]][current_95[[1]] == 0] <- NA

setwd(top_wd)
# pdf(file = "./outputs/figures/regrowth potential 2020 uncertainty.pdf",
#     width = 4.76 * 7.24 / 4.76,
#     height = 1.8 * 7.24 / 4.76)
png(file = "./outputs/figures/regrowth potential 2020 uncertainty.png",
    width = 4.76 * 7.24 / 4.76,
    height = 1.8 * 7.24 / 4.76,
    units = "in", res = 440)
library(colorspace)
ggplot() +
  geom_stars(data = current_95,
             downsample = 20,
             aes(x = x, y = y, fill = c30_2020_current_trop_for_and_sav_biomes_hi_fullmap.tif)) +
  scale_fill_gradientn(values = c(0, 1),
                       colours = c("grey95", "brown"),
                       na.value = "white") +
  guides(fill = guide_colorbar(title="Uncertainty in current regrowth\npotential (95% CI, Mg/ha/yr)",
                               barwidth = 12,
                               title.vjust = 1.7,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,-0.15),
        legend.direction = "horizontal") +
  theme(plot.margin=grid::unit(c(0,0,5,0), "mm")) +
  ylim(st_bbox(world)[2] * .7, st_bbox(world)[4] * .58) +
  xlim(st_bbox(world)[1] * .9, st_bbox(world)[3] * 1) +
  coord_equal() +
  theme_blank() +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA)
dev.off()

rm(current_95)
gc()


# Uncertainty: diff
setwd(raster_outputs_wd)
diff_lo <- read_stars("c30_2020_diff_trop_for_and_sav_biomes_lo_fullmap.tif",
                      proxy = F) + 0

diff_hi <- read_stars("c30_2020_diff_trop_for_and_sav_biomes_hi_fullmap.tif",
                      proxy = F) + 0

diff_95 <- diff_hi
diff_95[[1]] <- diff_hi[[1]] - diff_lo[[1]]
rm(diff_lo, diff_hi)
gc()

# Want to plot pixels without vales  as white not grey
diff_95[[1]][diff_95[[1]] == 0] <- NA


setwd(top_wd)
# pdf(file = "./outputs/figures/rewilding potential 2020 uncertainty.pdf",
#     width = 4.76 * 7.24 / 4.76,
#     height = 1.8 * 7.24 / 4.76)
png(file = "./outputs/figures/rewilding potential 2020 uncertainty.png",
    width = 4.76 * 7.24 / 4.76,
    height = 1.8 * 7.24 / 4.76,
    units = "in", res = 440)
library(colorspace)
ggplot() +
  geom_stars(data = diff_95,
             downsample = 20,
             aes(x = x, y = y, fill = c30_2020_diff_trop_for_and_sav_biomes_hi_fullmap.tif)) +
  scale_fill_gradientn(values = c(0, 1),
                       colours = c("grey95", "brown"),
                       na.value = "white") +
  guides(fill = guide_colorbar(title="Uncertainty in carbon accumulation potential of biodiversity\n& connectivity restoration (95% CI, Mg/ha/yr)",
                               barwidth = 12,
                               title.vjust = 1.7,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,-0.15),
        legend.direction = "horizontal") +
  theme(plot.margin=grid::unit(c(0,0,5,0), "mm")) +
  ylim(st_bbox(world)[2] * .7, st_bbox(world)[4] * .58) +
  xlim(st_bbox(world)[1] * .9, st_bbox(world)[3] * 1) +
  coord_equal() +
  theme_blank() +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA)
dev.off()

rm(diff_95)
gc()

















# Metrics ----------------------------------------------------------------------

# Calculating restoration potential in restoration sites ---------------

# For restoration sites, calculate potential carbon accumulation under current
# dispersal disruption and further carbon accumulation potential under seed
# dispersal disruption reversal

ncs_potential <- tibble(biome = c("trop_for", "trop_for", "trop_for_and_sav", "trop_for_and_sav"),
                        type = c("current", "diff", "current", "diff"),
                        mean = NA,
                        lo = NA,
                        hi = NA)

# First, consider tropical forest sites alone
setwd(raster_wd)
ncs_restore_sites <- read_stars("ncs_restore_sites.tif", proxy = F)

# Will mask to just tropical forest
ncs_restore_sites[[1]] <- ncs_restore_sites[[1]] * read_stars("trop_for_biomes.tif", proxy = F)[[1]]

# Get the number of total hectares
trop_for_restore_hectares <- ncs_restore_sites[[1]] %>% sum(na.rm = T) * proj_crs_cell_size


# trop_for current
ncs_ind <- which(ncs_potential$biome == "trop_for" & ncs_potential$type == "current")
setwd(top_wd)
ncs_potential$mean[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_current_trop_for_and_sav_biomes.tif",
                                           proxy = F)[[1]] *
                                  ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

ncs_potential$lo[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_current_trop_for_and_sav_biomes_lo_fullmap.tif",
                                         proxy = F)[[1]] *
                                ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

ncs_potential$hi[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_current_trop_for_and_sav_biomes_hi_fullmap.tif",
                                         proxy = F)[[1]] *
                                ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

# trop_for diff
ncs_ind <- which(ncs_potential$biome == "trop_for" & ncs_potential$type == "diff")
setwd(top_wd)
ncs_potential$mean[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_diff_trop_for_and_sav_biomes.tif",
                                           proxy = F)[[1]] *
                                  ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

ncs_potential$lo[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_diff_trop_for_and_sav_biomes_lo_fullmap.tif",
                                         proxy = F)[[1]] *
                                ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

ncs_potential$hi[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_diff_trop_for_and_sav_biomes_hi_fullmap.tif",
                                         proxy = F)[[1]] *
                                ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()



# Second, consider tropical forest and savanna sites
setwd(raster_wd)
ncs_restore_sites <- read_stars("ncs_restore_sites.tif", proxy = F)

# Will mask to tropical forest and savanna
ncs_restore_sites[[1]] <- ncs_restore_sites[[1]] * read_stars("trop_for_and_sav_biomes.tif", proxy = F)[[1]]

# Get the number of total hectares
trop_for_and_sav_restore_hectares <- ncs_restore_sites[[1]] %>% sum(na.rm = T) * proj_crs_cell_size


# trop_for current
ncs_ind <- which(ncs_potential$biome == "trop_for_and_sav" & ncs_potential$type == "current")
setwd(top_wd)
ncs_potential$mean[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_current_trop_for_and_sav_biomes.tif",
                                           proxy = F)[[1]] *
                                  ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

ncs_potential$lo[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_current_trop_for_and_sav_biomes_lo_fullmap.tif",
                                         proxy = F)[[1]] *
                                ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

ncs_potential$hi[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_current_trop_for_and_sav_biomes_hi_fullmap.tif",
                                         proxy = F)[[1]] *
                                ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

# trop_for diff
ncs_ind <- which(ncs_potential$biome == "trop_for_and_sav" & ncs_potential$type == "diff")
setwd(top_wd)
ncs_potential$mean[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_diff_trop_for_and_sav_biomes.tif",
                                           proxy = F)[[1]] *
                                  ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

ncs_potential$lo[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_diff_trop_for_and_sav_biomes_lo_fullmap.tif",
                                         proxy = F)[[1]] *
                                ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()

ncs_potential$hi[ncs_ind] <- (read_stars("./data/spatial data/raster_outputs/c30_2020_diff_trop_for_and_sav_biomes_hi_fullmap.tif",
                                         proxy = F)[[1]] *
                                ncs_restore_sites[[1]]) %>%
  sum(na.rm = T) * proj_crs_cell_size
gc()


# ncs_potential <- ncs_potential %>%
#   mutate(mean_per_ha = mean / ,
#          lo_per_ha = NA,
#          hi_per_ha = NA)

ncs_potential_pg <- ncs_potential %>%
  mutate(mean = mean / 1000000000,
         lo = lo / 1000000000,
         hi = hi / 1000000000)



setwd(top_wd)
png(file = "./outputs/figures/ncs summary.png",
    width = 180,
    height = 80,
    units = "mm", res = 440)
par(mar = c(4.1, 13, 4.1, 3.5))

yy1 <- 0.55
yy2 <- 2
yy3 <- 3.95

cex_units <- 0.9
cex_ncs <- 0.9

y_lab_diff <- 0.05
lwd_bar <- 20
lwd_ci <- lwd_bar / 10
col_ci <- "grey20"
  col_bar <- "grey70"


  plot(NA,
       ylim = c(4.75,-0.25),
       xlim = c(0, 2),
       las = 1,
       xlab = "",
       ylab = "",
       frame = F,
       xaxt = "n",
       yaxt = "n")
  axis(3, at = c(0, 0.5, 1, 1.5, 2),
       labels = c(0, 0.5, 1, 1.5, 2),
       mgp = c(3, 0.7, 0))

  text(0, -1.28, "PgC / yr", pos = 2, xpd = T, cex = cex_units)
  mtext("Climate mitigation potential in potential tropical restoration areas",
        side = 3,
        line = 2.5,
        font = 2)

  text("Current natural regrowth potential",
       x = 0,
       y = yy1 + y_lab_diff,
       pos = 2,
       xpd = T,
       cex = cex_ncs)

  segments(y0 = yy1,
           x0 = 0,
           x1 = ncs_potential_pg$mean[3],
           lwd = lwd_bar,
           lend = "butt",
           col = col_bar)

  arrows(y0 = yy1,
         x0 = ncs_potential_pg$lo[3],
         x1 = ncs_potential_pg$hi[3],
         lwd = lwd_ci,
         angle = 90,
         col = col_ci,
         code = 3,
         length = 0.04,
         lend = "butt")


  text("Regrowth potential if seed\ndispersal was not disrupted",
       x = 0,
       y = yy2 + y_lab_diff,
       pos = 2,
       xpd = T,
       cex = cex_ncs)

  segments(y0 = yy2,
           x0 = 0,
           x1 = ncs_potential_pg$mean[3] + ncs_potential_pg$mean[4],
           lwd = lwd_bar,
           lend = "butt",
           col = col_bar)

  arrows(y0 = yy2,
         x0 = ncs_potential_pg$lo[3] + ncs_potential_pg$lo[4],
         x1 = ncs_potential_pg$hi[3] + ncs_potential_pg$hi[4],
         lwd = lwd_ci,
         angle = 90,
         col = col_ci,
         code = 3,
         length = 0.04,
         lend = "butt",
         xpd = T)


  text("Carbon accumulation potential of \nbiodiversity and connectivity restoration",
       x = 0,
       y = yy3 + y_lab_diff,
       pos = 2,
       xpd = T,
       cex = cex_ncs)

  segments(y0 = yy3,
           x0 = 0,
           x1 = ncs_potential_pg$mean[4],
           lwd = lwd_bar,
           lend = "butt",
           col = col_bar)

  arrows(y0 = yy3,
         x0 = ncs_potential_pg$lo[4],
         x1 = ncs_potential_pg$hi[4],
         lwd = lwd_ci,
         angle = 90,
         col = col_ci,
         code = 3,
         length = 0.04,
         lend = "butt",
         xpd = T)

  to_Mg_per_ha <- 1000000000 / trop_for_and_sav_restore_hectares

  axis(1, at = c(0, 1.5, 3, 4.5) / to_Mg_per_ha,
       labels = c(0, 1.5, 3, 4.5),
       xpd = T,
       mgp = c(3, 0.5, 0),
       tck = -0.05)
  text(0, 5.7, "MgC / ha / yr", pos = 2, xpd = T, cex = cex_units)
  mtext("Average per-hectare climate mitigation potential",
        side = 1,
        line = 2.2,
        font = 2)

  dev.off()

# par(mar = c(1.1, 11.1, 5.1, 2.1))
#
# y_lab_diff <- 0.05
# lwd_bar <- 20
# lwd_ci <- lwd_bar / 8
# col_ci <- "grey20"
#   col_bar <- "grey70"
#
#   plot(NA,
#        ylim = c(4,0.5),
#        xlim = c(0, 2),
#        las = 1,
#        xlab = "",
#        ylab = "",
#        frame = F,
#        xaxt = "n",
#        yaxt = "n")
#   #axis(3, at = c(0, 0.6, 1.2))
#   axis(3, at = c(0, 0.5, 1, 1.5, 2))
#   mtext("Climate mitigation potential (PgC / yr)\nin potential tropical restoration areas (Griscom et al. 2017)",
#         side = 3,
#         line = 3,
#         font = 2)
#
#   text("Current natural\nregrowth potential",
#        x = 0,
#        y = 1 + y_lab_diff,
#        pos = 2,
#        xpd = T)
#
#   segments(y0 = 1,
#            x0 = 0,
#            x1 = ncs_potential_pg$mean[3],
#            lwd = lwd_bar,
#            lend = "butt",
#            col = col_bar)
#
#   arrows(y0 = 1,
#          x0 = ncs_potential_pg$lo[3],
#          x1 = ncs_potential_pg$hi[3],
#          lwd = lwd_ci,
#          angle = 90,
#          col = col_ci,
#          code = 3,
#          length = 0.04,
#          lend = "butt")
#
#
#   text("Regrowth potential if seed\ndispersal was not disrupted",
#        x = 0,
#        y = 2 + y_lab_diff,
#        pos = 2,
#        xpd = T)
#
#   segments(y0 = 2,
#            x0 = 0,
#            x1 = ncs_potential_pg$mean[3] + ncs_potential_pg$mean[4],
#            lwd = lwd_bar,
#            lend = "butt",
#            col = col_bar)
#
#   arrows(y0 = 2,
#          x0 = ncs_potential_pg$lo[3] + ncs_potential_pg$lo[4],
#          x1 = ncs_potential_pg$hi[3] + ncs_potential_pg$hi[4],
#          lwd = lwd_ci,
#          angle = 90,
#          col = col_ci,
#          code = 3,
#          length = 0.04,
#          lend = "butt",
#          xpd = T)
#
#
#   yval_3 <- 3.5
#   text("Carbon accumulation\npotential of biodiversity\nand connectivity restoration",
#        x = 0,
#        y = yval_3 + y_lab_diff,
#        pos = 2,
#        xpd = T)
#
#   segments(y0 = yval_3,
#            x0 = 0,
#            x1 = ncs_potential_pg$mean[4],
#            lwd = lwd_bar,
#            lend = "butt",
#            col = col_bar)
#
#   arrows(y0 = yval_3,
#          x0 = ncs_potential_pg$lo[4],
#          x1 = ncs_potential_pg$hi[4],
#          lwd = lwd_ci,
#          angle = 90,
#          col = col_ci,
#          code = 3,
#          length = 0.04,
#          lend = "butt",
#          xpd = T)
#
#
#
#   # Alternate on a per-hectare basis
#
#
#   ncs_potential_mg_ha <- ncs_potential %>%
#     mutate(mean = mean / c(trop_for_restore_hectares, trop_for_restore_hectares, trop_for_and_sav_restore_hectares, trop_for_and_sav_restore_hectares),
#            lo = lo / c(trop_for_restore_hectares, trop_for_restore_hectares, trop_for_and_sav_restore_hectares, trop_for_and_sav_restore_hectares),
#            hi = hi / c(trop_for_restore_hectares, trop_for_restore_hectares, trop_for_and_sav_restore_hectares, trop_for_and_sav_restore_hectares))
#
#
#   par(mar = c(1.1, 11.1, 5.1, 2.1))
#
#   y_lab_diff <- 0.05
#   lwd_bar <- 20
#   lwd_ci <- lwd_bar / 8
#   col_ci <- "grey20"
#     col_bar <- "grey70"
#
#
#     plot(NA,
#          ylim = c(4,0.5),
#          xlim = c(0, 4.5),
#          las = 1,
#          xlab = "",
#          ylab = "",
#          frame = F,
#          xaxt = "n",
#          yaxt = "n")
#     #axis(3, at = c(0, 0.6, 1.2))
#     axis(3, at = seq(0, 4.5, by = 1.5))
#     mtext("Average per-hectare climate mitigation potential (MgC / ha / yr)\nin potential tropical restoration areas (Griscom et al. 2017)",
#           side = 3,
#           line = 3,
#           font = 2)
#
#     text("Current natural\nregrowth potential",
#          x = 0,
#          y = 1 + y_lab_diff,
#          pos = 2,
#          xpd = T)
#
#     segments(y0 = 1,
#              x0 = 0,
#              x1 = ncs_potential_mg_ha$mean[3],
#              lwd = lwd_bar,
#              lend = "butt",
#              col = col_bar)
#
#     arrows(y0 = 1,
#            x0 = ncs_potential_mg_ha$lo[3],
#            x1 = ncs_potential_mg_ha$hi[3],
#            lwd = lwd_ci,
#            angle = 90,
#            col = col_ci,
#            code = 3,
#            length = 0.04,
#            lend = "butt")
#
#
#     text("Regrowth potential if seed\ndispersal was not disrupted",
#          x = 0,
#          y = 2 + y_lab_diff,
#          pos = 2,
#          xpd = T)
#
#     segments(y0 = 2,
#              x0 = 0,
#              x1 = ncs_potential_mg_ha$mean[3] + ncs_potential_mg_ha$mean[4],
#              lwd = lwd_bar,
#              lend = "butt",
#              col = col_bar)
#
#     arrows(y0 = 2,
#            x0 = ncs_potential_mg_ha$lo[3] + ncs_potential_mg_ha$lo[4],
#            x1 = ncs_potential_mg_ha$hi[3] + ncs_potential_mg_ha$hi[4],
#            lwd = lwd_ci,
#            angle = 90,
#            col = col_ci,
#            code = 3,
#            length = 0.04,
#            lend = "butt",
#            xpd = T)
#
#
#     yval_3 <- 3.5
#     text("Carbon accumulation\npotential of biodiversity\nand connectivity restoration",
#          x = 0,
#          y = yval_3 + y_lab_diff,
#          pos = 2,
#          xpd = T)
#
#     segments(y0 = yval_3,
#              x0 = 0,
#              x1 = ncs_potential_mg_ha$mean[4],
#              lwd = lwd_bar,
#              lend = "butt",
#              col = col_bar)
#
#     arrows(y0 = yval_3,
#            x0 = ncs_potential_mg_ha$lo[4],
#            x1 = ncs_potential_mg_ha$hi[4],
#            lwd = lwd_ci,
#            angle = 90,
#            col = col_ci,
#            code = 3,
#            length = 0.04,
#            lend = "butt",
#            xpd = T)
#



# Another way to explore spatial variation at the country/region level and
# associated uncertainty surrounding the lost regrowth potential or the
# restoration potential.

# library("tidyverse")
# library(stars)
# library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Load country polygons
countries <- ne_countries(scale = "medium", returnclass = "sf")


# Load in the ncs sites
# Consider tropical forest and savanna sites
setwd(raster_wd)
ncs_restore_sites <- read_stars("ncs_restore_sites.tif", proxy = F)

# Will mask to tropical forest and savanna
ncs_restore_sites[[1]] <- ncs_restore_sites[[1]] * read_stars("trop_for_and_sav_biomes.tif", proxy = F)[[1]]

# Want to convert to NA all the pixels that are not potential restoration sites
ncs_restore_sites[[1]][ncs_restore_sites[[1]] == 0] <- NA



# Set up for adding these stats to the countries tibble
countries <- countries %>%
  st_transform(st_crs(diff_ncs_sites)) %>%
  mutate(sum_diff = NA,
         mean_diff = NA,
         sum_diff_lo = NA,
         mean_diff_lo = NA,
         sum_diff_hi = NA,
         mean_diff_hi = NA)

gc()




# Pull in the diff 2020 raster and turn non-restoration sites to NA
setwd(raster_outputs_wd)
diff_ncs_sites <- read_stars("c30_2020_diff_trop_for_and_sav_biomes.tif",
                             proxy = F)
diff_ncs_sites[[1]] <- diff_ncs_sites[[1]] * ncs_restore_sites[[1]]
rm(ncs_restore_sites)
gc()

# This was really slow using the stars package, so will try
# the exactextractr package # install.packages("exactextractr")
library("exactextractr")
library("raster")
st_as_raster <- function(rstars){
  rext <- st_bbox(rstars)
  raster(t(rstars[[1]]), xmn = rext[1], xmx = rext[3],
         ymn = rext[2], ymx=rext[4],
         crs = st_crs(rstars)$proj4string)
}
diff_ncs_sites_raster <- st_as_raster(diff_ncs_sites)


countries$mean_diff <- exact_extract(diff_ncs_sites_raster,
                                     countries %>% dplyr::select(name),
                                     fun = "mean")
gc()

proj_crs_cell_size <- 708264.2 / 10000
countries$sum_diff <- exact_extract(diff_ncs_sites_raster,
                                    countries %>% dplyr::select(name),
                                    fun = "sum") * proj_crs_cell_size
gc()

rm(diff_ncs_sites, diff_ncs_sites_raster)




# Pull in the diff_lo 2020 raster and turn non-restoration sites to NA
setwd(raster_outputs_wd)
diff_lo_ncs_sites <- read_stars("c30_2020_diff_trop_for_and_sav_biomes_lo_fullmap.tif",
                                proxy = F)
diff_lo_ncs_sites[[1]] <- diff_lo_ncs_sites[[1]] * ncs_restore_sites[[1]]

gc()

diff_lo_ncs_sites_raster <- st_as_raster(diff_lo_ncs_sites)

gc()
countries$mean_diff_lo <- exact_extract(diff_lo_ncs_sites_raster,
                                        countries %>% dplyr::select(name),
                                        fun = "mean")
gc()

countries$sum_diff_lo <- exact_extract(diff_lo_ncs_sites_raster,
                                       countries %>% dplyr::select(name),
                                       fun = "sum") * proj_crs_cell_size
gc()

rm(diff_lo_ncs_sites, diff_lo_ncs_sites_raster)



# Pull in the diff_hi 2020 raster and turn non-restoration sites to NA
setwd(raster_outputs_wd)
diff_hi_ncs_sites <- read_stars("c30_2020_diff_trop_for_and_sav_biomes_hi_fullmap.tif",
                                proxy = F)
diff_hi_ncs_sites[[1]] <- diff_hi_ncs_sites[[1]] * ncs_restore_sites[[1]]

gc()

diff_hi_ncs_sites_raster <- st_as_raster(diff_hi_ncs_sites)

gc()
countries$mean_diff_hi <- exact_extract(diff_hi_ncs_sites_raster,
                                        countries %>% dplyr::select(name),
                                        fun = "mean")
gc()

countries$sum_diff_hi <- exact_extract(diff_hi_ncs_sites_raster,
                                       countries %>% dplyr::select(name),
                                       fun = "sum") * proj_crs_cell_size
gc()

rm(diff_hi_ncs_sites, diff_hi_ncs_sites_raster)











# Write out a table
countries %>%
  st_drop_geometry() %>%
  dplyr::select(name, sum_diff,
                sum_diff_lo,
                sum_diff_hi,
                mean_diff,
                mean_diff_lo,
                mean_diff_hi) %>%
  mutate(sum_diff = round(sum_diff),
         mean_diff = round(mean_diff, 2),
         sum_diff_lo = round(sum_diff_lo),
         mean_diff_lo = round(mean_diff_lo, 2),
         sum_diff_hi = round(sum_diff_hi),
         mean_diff_hi = round(mean_diff_hi, 2)) %>%
  dplyr::filter(sum_diff > 0) %>%
  dplyr::arrange(desc(sum_diff)) %>%
  write.csv(file = "~/Downloads/country_summary.csv")



nd_diff <- read_stars("~/Downloads/c30_2020_diff_nd_trop_for_and_sav_biomes.tif")
nd <- read_stars("~/Downloads/c30_2020_nd_trop_for_and_sav_biomes.tif")


setwd(top_wd)
# pdf(file = "./outputs/figures/nd regrowth potential 2020.pdf",
#     width = 4.76 * 7.24 / 4.76,
#     height = 1.8 * 7.24 / 4.76)
png(file = "./outputs/figures/nd regrowth potential 2020.png",
    width = 4.76 * 7.24 / 4.76,
    height = 1.8 * 7.24 / 4.76,
    units = "in", res = 440)

library(colorspace)
ggplot() +
  geom_stars(data = nd,
             downsample = 20,
             aes(x = x, y = y, fill = c30_2020_nd_trop_for_and_sav_biomes.tif)) +
  scale_fill_gradientn(values = c(0, 0.4, 1),
                       colours = c("lightyellow", "forestgreen", "darkgreen"),
                       na.value = "white") +
  guides(fill = guide_colorbar(title="Regrowth potential from\nabiotic-only model (Mg/ha/yr)",
                               barwidth = 12,
                               title.vjust = 1.7,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,-0.15),
        legend.direction = "horizontal") +
  theme(plot.margin=grid::unit(c(0,0,5,0), "mm")) +
  ylim(st_bbox(world)[2] * .7, st_bbox(world)[4] * .58) +
  xlim(st_bbox(world)[1] * .9, st_bbox(world)[3] * 1) +
  coord_equal() +
  theme_blank() +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA)
dev.off()





png(file = "./outputs/figures/regrowth potential difference from abiotic only.png",
    width = 4.76 * 7.24 / 4.76,
    height = 1.8 * 7.24 / 4.76,
    units = "in", res = 440)

library(colorspace)
ggplot() +
  geom_stars(data = nd_diff,
             downsample = 20,
             aes(x = x, y = y, fill = c30_2020_diff_nd_trop_for_and_sav_biomes.tif)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,
                                   na.value = "white",
                                   p1 = 0.6,
                                   p2 = 0.6,
                                   p3 = 0.6,
                                   p4 = 0.6) +
  guides(fill = guide_colorbar(title="Difference between dispersal disruption-\ndependent and abiotic-only model (Mg/ha/yr)",
                               barwidth = 12,
                               title.vjust = 1.7,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(#legend.position="bottom",
    legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,-0.15),
        legend.direction = "horizontal") +
  theme(plot.margin=grid::unit(c(0,0,5,0), "mm")) +
  ylim(st_bbox(world)[2] * .7, st_bbox(world)[4] * .58) +
  xlim(st_bbox(world)[1] * .9, st_bbox(world)[3] * 1) +
  coord_equal() +
  theme_blank() +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA)
dev.off()

