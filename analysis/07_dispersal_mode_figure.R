source(list.files("./R", full.names = T))

ipak(c("tidyverse", "RCurl"))


# Load in the sPlot data -------------------------------------------------------
# https://idata.idiv.de/ddm/Data/ShowData/3474?version=54
list.files("./data") # This has been shared as an RData objected
load("./data/sPlotOpen/sPlotOpen.RData") #


# Munging the plot data --------------------------------------------------------

# Rename
splot_dat <- DT2.oa %>%
  left_join(header.oa)
rm(DT2.oa)


# Subset to forest only
splot_dat$Biome %>% table()
splot_dat <- splot_dat %>%
  filter(Forest,
         is_forest)

# Clean up names

# Clean up names
# Want to clean this up with TNRS
splot_sp <- tibble(sp = unique(splot_dat$Species),
                   gen = word(unique(splot_dat$Species), 1))


# This takes a couple minutes
# Takes a while, so will save
if(!"splot_tnrs.csv" %in% list.files("./data")){
  splot_tnrs <- TNRS::TNRS(splot_sp$sp)
  write.csv(splot_tnrs, "./data/splot_tnrs.csv")
} else{
  splot_tnrs <- read.csv("./data/splot_tnrs.csv", row.names = 1) %>% tibble()
}

# Put this on the splot_sp dataframe
splot_sp <- splot_sp %>%
  left_join(splot_tnrs, by = c("sp" = "Name_submitted")) %>%
  mutate(tnrs_sp = Accepted_species,
         tnrs_fam = Accepted_family) %>%
  dplyr::select("sp",
                "gen",
                "tnrs_sp",
                "tnrs_fam") %>%
  mutate(tnrs_sp = ifelse(tnrs_sp == "",
                          NA,
                          tnrs_sp))

# Then join this back to splot_dat
splot_dat <- splot_dat %>%
  left_join(splot_sp, by = c("Species" = "sp"))



# Assign dispersal mode --------------------------------------------------------

# Get dispersal mode info
mode_dat <- getURL("https://raw.githubusercontent.com/evancf/classifying-dispersal-mode/main/data/mode_dat.csv")
mode_dat <- read.csv(text = mode_dat, row.names = 1) %>%
  tibble()

# Want to clean this up with TNRS
mode_sp <- tibble(sp = unique(mode_dat$sp),
                  gen = word(unique(mode_dat$sp), 1))

# This takes a couple minutes
if(!"mode_tnrs.csv" %in% list.files("./data")){
  mode_tnrs <- TNRS::TNRS(mode_sp$sp)
  write.csv(mode_tnrs, "./data/mode_tnrs.csv")
} else{
  mode_tnrs <- read.csv("./data/mode_tnrs.csv", row.names = 1) %>% tibble()
}

# Put this on the mode_sp dataframe
mode_sp <- mode_sp %>%
  mutate(tnrs_sp = mode_tnrs$Accepted_species,
         tnrs_fam = mode_tnrs$Accepted_family) %>%
  mutate(tnrs_sp = ifelse(tnrs_sp == "",
                          NA,
                          tnrs_sp))


# Then join this back to mode_dat
mode_dat <- mode_dat %>%
  left_join(mode_sp)




# See the coverage in splot

(unique(splot_dat$tnrs_sp) %in% mode_dat$tnrs_sp) %>% mean()
((splot_dat %>%
    filter(Latitude > -23.5,
           Latitude < 23.5) %>%
    pull(tnrs_sp)) %in% mode_dat$tnrs_sp) %>% mean()

colnames(splot_dat)



# Try to merge the plot and mode data by species

mode_mean_sp <- mode_dat %>%
  group_by(tnrs_sp) %>%
  summarise(bio_sp = mean(biotic))

mode_mean_gen <- mode_mean_sp %>%
  ungroup() %>%
  mutate(tnrs_gen = word(tnrs_sp, 1)) %>%
  group_by(tnrs_gen) %>%
  summarise(bio_gen = mean(bio_sp))

mode_mean_fam <- mode_mean_gen %>%
  ungroup() %>%
  left_join(unique(mode_dat %>% dplyr::select(gen, tnrs_fam)),
            by = c("tnrs_gen" = "gen")) %>%
  group_by(tnrs_fam) %>%
  summarise(bio_fam = mean(bio_gen))


hist(mode_mean_sp$bio_sp)
hist(mode_mean_gen$bio_gen)
hist(mode_mean_fam$bio_fam)

# Join these together

splot_sp_with_mode <- splot_sp %>%
  rename(tnrs_gen = gen) %>%
  left_join(mode_mean_sp) %>%
  left_join(mode_mean_gen) %>%
  left_join(mode_mean_fam) %>%
  mutate(biotic = bio_sp) %>%
  mutate(biotic = ifelse(is.na(biotic),
                         bio_gen,
                         biotic)) %>%
  mutate(biotic = ifelse(is.na(biotic),
                         bio_fam,
                         biotic))

# Now join to the splot_dat
splot_dat <- splot_dat %>%
  left_join(dplyr::select(splot_sp_with_mode, tnrs_sp, biotic) %>%
              filter(!is.na(tnrs_sp),
                     !duplicated(tnrs_sp)))

mean(splot_dat$biotic, na.rm=T)


# Get the plot level averages

splot_mean <- splot_dat %>%
  group_by(PlotObservationID) %>%
  summarize(biotic = mean(biotic, na.rm = T),
            Latitude = first(Latitude),
            Longitude = first(Longitude),
            Naturalness = first(Naturalness),
            Grassland = first(Grassland),
            Shrubland = first(Shrubland),
            Wetland = first(Wetland),
            Cover_tree_layer = first(Cover_tree_layer),
            Dataset = first(Dataset))

splot_mean_forest <- splot_mean %>%
  filter(Naturalness == "Natural",
         !Grassland,
         !Shrubland,
         !Wetland)




# Also want to include BIEN data -----------------------------------------------

# Because these include info on dbh, we can
ipak("BIEN")

# Pull data from BIEN
if(!"stem_dat.csv" %in% list.files("./data")){
  stem1 <- BIEN_stem_sampling_protocol("0.1 ha  transect, stems >= 2.5 cm dbh")
  stem2 <- BIEN_stem_sampling_protocol("0.1 ha  transect, stems >= 2.5 cm dbh (Warning: species + count of individuals + stem DBHs recorded once per subplot. Individuals not linked to stems)")

  stem_dat <- bind_rows(stem1, stem2) %>%
    tibble() %>%
    mutate(tnrs_gen = word(scrubbed_species_binomial, 1))
  stem_dat <- stem_dat %>%
    left_join(TNRS::TNRS(stem_dat$tnrs_gen %>% unique()) %>% dplyr::select(Name_submitted, Accepted_family),
              by = c("tnrs_gen" = "Name_submitted")) %>%
    mutate(tnrs_fam = Accepted_family) %>%
    dplyr::select(-Accepted_family)
  write.csv(stem_dat, "./data/stem_dat.csv")
} else{
  stem_dat <- read.csv("./data/stem_dat.csv", row.names = 1) %>% tibble()
}


# Get ssd data from BIEN as well

if(!"ssd_dat.csv" %in% list.files("./data")){
  ssd_dat <- BIEN_trait_trait("stem wood density") %>% tibble()
  write.csv(ssd_dat, "./data/ssd_dat.csv")
} else{
  ssd_dat <- read.csv("./data/ssd_dat.csv", row.names = 1) %>% tibble()
}

# Not sure if plant growth form would  be useful...
# pgf_dat <- BIEN_trait_trait("whole plant growth form") %>% tibble()

ssd_mean_sp <- ssd_dat %>%
  group_by(scrubbed_species_binomial) %>%
  summarise(ssd_sp = mean(as.numeric(trait_value, na.rm = T))) %>%
  rename(tnrs_sp = scrubbed_species_binomial)

ssd_mean_gen <- ssd_mean_sp %>%
  ungroup() %>%
  mutate(tnrs_gen = word(tnrs_sp, 1)) %>%
  group_by(tnrs_gen) %>%
  summarise(ssd_gen = mean(ssd_sp, na.rm = T))

ssd_mean_fam <- ssd_mean_gen %>%
  ungroup() %>%
  left_join(TNRS::TNRS(ssd_mean_gen$tnrs_gen %>% unique()) %>% dplyr::select(Name_submitted, Accepted_family),
            by = c("tnrs_gen" = "Name_submitted")) %>%
  mutate(tnrs_fam = Accepted_family) %>%
  group_by(tnrs_fam) %>%
  summarise(ssd_fam = mean(ssd_gen, na.rm = T))


# Some renaming
stem_dat <- stem_dat %>%
  mutate(tnrs_sp = scrubbed_species_binomial,
         tnrs_gen = word(scrubbed_species_binomial, 1))

(stem_dat$tnrs_sp %in% ssd_mean_sp$tnrs_sp |
    stem_dat$tnrs_gen %in% ssd_mean_gen$tnrs_gen |
    stem_dat$tnrs_fam %in% ssd_mean_fam$tnrs_fam) %>% table()

# Assign ssd to each species in the stem data

stem_sp <- stem_dat %>%
  filter(!duplicated(tnrs_sp))
stem_sp %>% dim()


stem_sp_with_ssd <- stem_sp %>%
  left_join(ssd_mean_sp) %>%
  left_join(ssd_mean_gen) %>%
  left_join(ssd_mean_fam) %>%
  mutate(mean_ssd = ssd_sp) %>%
  mutate(mean_ssd = ifelse(is.na(mean_ssd),
                           ssd_gen,
                           mean_ssd)) %>%
  mutate(mean_ssd = ifelse(is.na(mean_ssd),
                           ssd_fam,
                           mean_ssd))
# Look for ones that don't have a value
stem_sp_with_ssd %>%
  dplyr::select(scrubbed_species_binomial,
                ssd_sp, ssd_gen, ssd_fam) %>%
  pull(ssd_fam) %>% is.na() %>% table()


# Now join to the stem_dat
stem_dat <- stem_dat %>%
  left_join(dplyr::select(stem_sp_with_ssd, tnrs_sp, mean_ssd) %>%
              filter(!is.na(tnrs_sp),
                     !duplicated(tnrs_sp)))

average_ssd <- mean(stem_dat$mean_ssd, na.rm=T)


# Get biomass data for the BIEN plots ------------------------------------------

ipak("BIOMASS")

# A few more manipulations
short_stem_dat <- stem_dat %>%
  mutate(mean_ssd = ifelse(is.na(mean_ssd),
                           average_ssd,
                           mean_ssd)) %>%  # If we don't have SSD data, we will assume overall mean values
  filter(!plot_name %in% c("PROVIDEN")) # There's an island plots that won't work with a height model raster

# bad_indices <- FILL IN INDICES IF THERE ARE ERRORS IN THE computeAGB function
# short_stem_dat$plot_name[bad_indices]

short_stem_dat <- short_stem_dat %>%
  mutate(agb = computeAGB(D = short_stem_dat$stem_dbh_cm,
                          WD = short_stem_dat$mean_ssd,
                          coord = short_stem_dat[, c("longitude", "latitude")])) %>%
  mutate(agb = agb * individual_count)



# Also want to get dispersal mode information added to this

short_stem_dat <- short_stem_dat %>%
  left_join(mode_mean_sp) %>%
  left_join(mode_mean_gen) %>%
  left_join(mode_mean_fam) %>%
  mutate(biotic = bio_sp) %>%
  mutate(biotic = ifelse(is.na(biotic),
                         bio_gen,
                         biotic)) %>%
  mutate(biotic = ifelse(is.na(biotic),
                         bio_fam,
                         biotic))



# Now summarize to get an equivalent of splot_mean

bien_mean <- short_stem_dat %>%
  #ungroup() %>%
  group_by(plot_name) %>%
  summarize(sum_agb = sum(agb, na.rm = T),
            weighted_biotic_sum_agb = sum(biotic * agb, na.rm = T),
            plot_area_ha = dplyr::first(plot_area_ha),
            latitude = first(latitude),
            longitude = first(longitude)) %>%
  mutate(biotic = weighted_biotic_sum_agb / sum_agb)

# Find shared plots that are shared between the two datasets.
# Unfortunately, we don't have shared plot names
(splot_mean$PlotObservationID %in% bien_mean$plot_name) %>% table()

# But we do have the ability to match by long and lat
paste(splot_mean$Latitude, splot_mean$Longitude) %in% paste(bien_mean$latitude, bien_mean$longitude) %>% table()

splot_mean <- splot_mean %>%
  mutate(coord_id = paste(Latitude, Longitude, sep = "_"))

bien_mean <- bien_mean %>%
  mutate(coord_id = paste(latitude, longitude, sep = "_"))

shared_coords <- intersect(splot_mean$coord_id, bien_mean$coord_id)

splot_mean$coord_id %in% bien_mean$coord_id %>% table()
bien_mean$coord_id %in% splot_mean$coord_id %>% table()

# Get the subset of each data set that includes those coordinates
splot_set <- splot_mean %>%
  filter(coord_id %in% shared_coords) %>%
  rename(biotic_splot = biotic)
bien_set <- bien_mean %>%
  filter(coord_id %in% shared_coords) %>%
  rename(biotic_bien = biotic)
# Now join these together
shared_set <- inner_join(splot_set,
                         bien_set)

plot(biotic_splot ~ biotic_bien,
     data = shared_set)


mean(shared_set$biotic_bien, na.rm = T)
plot(biotic ~ latitude, data = bien_mean)



# Merge the bien data and splot data, recognizing that we need to remove the
# problem associated with the Salvias data present in sPlot

bien_mean %>% glimpse()
splot_mean_forest %>% glimpse()

splot_mean_to_use <- splot_mean_forest %>%
  filter(Dataset != "Salvias") %>%  # Remove the Salvias data
  dplyr::rename(plot_id = PlotObservationID,
                latitude = Latitude,
                longitude = Longitude) %>%
  dplyr::select(plot_id,
                latitude,
                longitude,
                biotic) %>%
  mutate_at(1, factor)

bien_mean_to_use <- bien_mean %>%
  dplyr::rename(plot_id = plot_name) %>% #  These are all Salvias sites...
  mutate(biotic = weighted_biotic_sum_agb / sum_agb) %>%
  dplyr::select(plot_id,
                latitude,
                longitude,
                biotic)


both_mean <- bind_rows(splot_mean_to_use,
                       bien_mean_to_use)

# # Write this out for later use
# write.csv(both_mean, "./data/both_mean.csv")

# Random forest-based map --------------------------------------------------

both_mean <- read.csv("./data/both_mean.csv", header = T, row.names = 1) %>%
  tibble()

# ipak(c("dismo", "maptools",
#        "geodata", "raster"))
ipak(c("raster"))


# # To get worldclim data on your local computer
# wc_dat <- worldclim_global(var = "bio", res = 5, path = tempfile())
#
# wd <- getwd()
# setwd("./data/worldclim")
# for(i in 1:length(names(wc_dat))){
#   raster::writeRaster(wc_dat[[i]],
#               filename = paste(names(wc_dat[[i]]),"tif", sep = "."))
# }
# setwd(wd)

wc_dat <- stack(list.files(path = "./data/worldclim", pattern = '\\.tif$', full.names = T))

bio_vals <- raster::extract(wc_dat,
                            both_mean[,c("longitude", "latitude")] %>%
                              dplyr::rename(lon = longitude,
                                            lat = latitude))

both_mean <- both_mean %>%
  bind_cols(bio_vals)

# Not sure what the dplyr approach to this is...
colnames(both_mean) <- gsub("wc2.1_5m_", "", colnames(both_mean))


ipak("tidymodels")
ipak("sf")
set.seed(4)
both_mean_sf <- both_mean %>%
  dplyr::select(-plot_id) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326)
both_split <- initial_split(both_mean_sf)
both_train <- training(both_split)
both_test <- testing(both_split)

both_rec <- recipe(biotic ~ .,
                   data = both_train %>% st_drop_geometry())

# Let's make sure this looks good
both_prep <- prep(both_rec)
juiced <- juice(both_prep)

tune_spec <- rand_forest(
  mtry = tune(),
  trees = tune(),
  min_n = tune()
) %>%
  set_mode("regression") %>%
  set_engine("ranger")


tune_wf <- workflow() %>%
  add_recipe(both_rec) %>%
  add_model(tune_spec)


# Spatial clustering cv
ipak("spatialsample")
set.seed(5)
both_folds <- spatial_clustering_cv(both_train, v = 5)
autoplot(both_folds)

# Here are some tuning steps. This takes a while
library("tictoc")
tic()
set.seed(6)
tune_res <- tune_grid(
  tune_wf,
  resamples = both_folds,
  grid = 50 # Should take about an hour
)
toc()

tune_res


tune_res %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  dplyr::select(mean, min_n, trees, mtry) %>%
  pivot_longer(min_n:mtry,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "rsq")

# Because the above takes a while, will just hard code in the values that work best
# mtry = 1, trees = 821, min_n = 5 .config = Preprocessor1_Model23
rf_grid <- grid_regular(
  mtry(range = c(1, 1)),
  trees(range = c(821, 821)),
  min_n(range = c(5, 5)),
  levels = 1
)

set.seed(7)
regular_res <- tune_grid(
  tune_wf,
  resamples = both_folds,
  grid = rf_grid
)

best_rsq <- select_best(regular_res, "rsq")

final_rf <- finalize_model(tune_spec,
                           best_rsq)

final_wf <- workflow() %>%
  add_recipe(both_rec) %>%
  add_model(final_rf)

final_res <- final_wf %>%
  last_fit(both_split)

# Get the R squared for this model on the test set
final_res %>%
  collect_metrics()
# Checking this makes sense
lm(both_test$biotic ~ (final_res %>%
                         collect_predictions() %>%
                         pull(.pred))) %>%
  summary()

set.seed(4)
for_pred <- finalize_workflow(tune_wf,
                              best_rsq) %>%
  fit(both_mean)



# Now let's predict on rasters
wc_newdata <- wc_dat
names(wc_newdata) <- gsub("wc2.1_5m_", "", names(wc_newdata))


# Tips from here: https://stackoverflow.com/questions/66787014/using-parsnip-model-to-predict-a-raster-in-r
get_pred_fun<-function(...){
  p<-predict(...)
  return(as.matrix(as.numeric(p[, 1, drop=T])))
}
biotic_map <- raster::predict(wc_newdata,
                              for_pred,
                              fun = get_pred_fun)

writeRaster(biotic_map, "./data/biotic_map.tif", overwrite=T)


# Lets now assess uncertainty --------------------------------------------------

# We will use bootstrapping to map model uncertainty
set.seed(4)
for(i in 1:10){
  for_pred_boot <- finalize_workflow(tune_wf,
                                     best_rsq) %>%
    fit(both_mean[sample(1:nrow(both_mean), replace = T),])


  biotic_map_boot <- raster::predict(wc_newdata,
                                     for_pred_boot,
                                     fun = get_pred_fun)

  #plot(biotic_map_boot)

  writeRaster(biotic_map_boot, paste0("./data/biotic_map_boot", i, ".tif"), overwrite=T)
}

# We will load these back in as a raster stack

boot_stack <- stack(list.files("./data/", pattern = "biotic_map_boot",
                               full.names = T))

boot_mean <- calc(boot_stack, fun = mean, na.rm = T)
boot_sd <- calc(boot_stack, fun = sd, na.rm = T)

boot_coef_var <- boot_sd/boot_mean

plot(boot_coef_var)

writeRaster(boot_coef_var, "./data/boot_coef_var.tif", overwrite=T)


# Now for plotting -----------------------------------------------------------

library(raster)
biotic_map <- raster("./data/biotic_map.tif", proxy = F)
boot_coef_var <- raster("./data/boot_coef_var.tif", proxy = F)
# Transform to match other maps
proj_crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
biotic_map <- projectRaster(biotic_map, crs = CRS(proj_crs))
boot_coef_var <- projectRaster(boot_coef_var, crs = CRS(proj_crs))

for_and_sav_biomes <- raster("./data/spatial data/raster/for_and_sav_biomes.tif", proxy = F) %>%
  projectRaster(to = biotic_map)

for_and_sav_biomes <- round(for_and_sav_biomes)
plot(biotic_map * for_and_sav_biomes)

# First get a continents polygon
# Want to match this stars raster
library("rnaturalearth")
sf_use_s2(FALSE)
library("stars")
world <- ne_countries() %>%
  st_as_sf() %>%
  filter(continent != "Antarctica") %>%
  st_union() %>%
  st_transform(crs = st_crs(for_and_sav_biomes %>% st_as_stars()))
focal_bbox <- st_bbox(for_and_sav_biomes %>% st_as_stars())
sf_use_s2(T)

# Also want to visualize the plot data as a scatter plot

# Deal with values equal to 0 and 1 based on maximum and minimum values
biotic_range <- both_mean$biotic[!both_mean$biotic %in% c(0, 1)] %>% range()

both_mean$biotic_logit <- ifelse(both_mean$biotic > (1-biotic_range[1]),
                                 (1-biotic_range[1]), both_mean$biotic)
both_mean$biotic_logit <- ifelse(both_mean$biotic_logit < biotic_range[1],
                                 biotic_range[1], both_mean$biotic_logit)
both_mean$biotic_logit <-   both_mean$biotic_logit %>%
  gtools::logit()

# Will use a gam to fit a curve
ipak("sjPlot")
gam_mod <- mgcv::gam(biotic_logit ~ s(latitude, k = 5),
                     data = both_mean)
gam_preds <- get_model_data(gam_mod, type = "pred",
                            transform = "plogis")


#pdf(
png(file = "./outputs/figures/biotic dispersal in forests and savannah.png",
    width = 4.76,
    height = 2.25,
    units = "in", res = 440)

plot_biotic_map <- biotic_map * for_and_sav_biomes

plot_biotic_df <- as.data.frame(plot_biotic_map, xy=T) %>%
  filter(!is.na(layer)) %>%
  mutate(layer = ifelse(layer == 0, NA, layer))


p1 <- ggplot() +
  geom_raster(data = plot_biotic_df, aes(fill = layer,
                                         x=x, y=y)) +
  scale_fill_viridis_c(limits = c(0, 1),
                       values = c(0,0.3, 0.5, 0.6, 0.7, 1),
                       direction = -1,
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0%", "25%", "50%", "75%", "100%"),) +
  coord_quickmap() +
  theme_void() +
  theme(plot.margin = margin(0.2,0,0,0, "cm")) +
  #ylim(-60, 90) +
  guides(fill = guide_colorbar(title="Animal seed dispersal prevalence",
                               barwidth = 10,
                               title.vjust = 1.14,
                               label.vjust = 1.6,
                               barheight = 0.4)) +
  theme(#legend.position="bottom",
    legend.title=element_text(size=10)) +
  theme(legend.position = c(.48,1.04),
        legend.direction = "horizontal") +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA)

# guides(fill = guide_colorbar(title="Biotic dispersal\nprevalence",
#                              barwidth = 12,
#                              title.vjust = 1.5,
#                              label.vjust = 1.6,
#                              barheight = 0.46)) +
#   theme(#legend.position="bottom",
#     legend.title=element_text(size=10)) +
#   theme(legend.position = c(.4,.925),
#         legend.direction = "horizontal")

my_inset <- function(){
  par(mar = c(0,0,0,0),
      mgp = c(1,0,0)
  )
  plot(NA,
       ylim = c(-60, 80),
       xlim = c(0, 1),
       frame = F,
       las = 1,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       cex = 0.5)
  axis(2, at = c(-60,0,60), labels = c("60°S  ", "0  ", "60°N  "), las = 1, tck = -0.005,
       cex.axis = 0.7)
  axis(1, at = c(0,0.5,1), labels = c("0%", "50%", "100%"), las = 1, tck = -0.005,
       cex.axis = 0.7)


  set.seed(4)
  points(jitter(latitude, 40) ~ biotic,
         data = both_mean,
         pch = 16,
         cex = 0.05,
         col = "grey60")#rgb(0,0,0,0.1))

  lines(gam_preds$latitude$x ~ plogis(gam_preds$latitude$predicted),
        col = "darkgreen")

  polygon(x = plogis(c(gam_preds$latitude$conf.low,
                       rev(gam_preds$latitude$conf.high))),
          y = c(gam_preds$latitude$x,
                rev(gam_preds$latitude$x)),
          col = rgb(34, 139, 34, maxColorValue = 255, alpha = 150), #2, 48, 32
          border = NA)
}


ipak("grid")
ipak("ggplotify")

p2 <- as.grob(~my_inset())

#grid.newpage()
grid.draw(p1)
vp = viewport(x=.18, y=.34, width=.135, height=0.3)#.275)
pushViewport(vp)
grid.draw(p2)
upViewport()


dev.off()





# Now for the uncertainty one
#pdf(
png(file = "./outputs/figures/biotic dispersal in forests and savannah uncertainty.png",
    width = 6,
    height = 2.5,
    units = "in", res = 440)

plot_boot_coef_var <- boot_coef_var * for_and_sav_biomes

plot_biotic_df <- as.data.frame(plot_boot_coef_var, xy=T) %>%
  filter(!is.na(layer)) %>%
  mutate(layer = ifelse(layer == 0, NA, layer))


ggplot() +
  geom_raster(data = plot_biotic_df, aes(fill = layer,
                                         x=x, y=y)) +
  scale_fill_viridis_c(limits = c(0, 0.2),
                       direction = -1) +
  coord_quickmap() +
  theme_void() +
  #ylim(-60, 90) +
  guides(fill = guide_colorbar(title="Biotic\ndispersal\nprevalence\nuncertainty\n(coef. var.)")) +
  geom_sf(data = world,
          col = "grey",
          linewidth = 0.1,
          fill = NA)

dev.off()



# # Alternative with color scale at bottom
#
# pdf(file = "./outputs/figures/biotic dispersal in forests and savannah.pdf",
#     width = 4.76,
#     height = 2.5)
#
# plot_biotic_map <- biotic_map * for_and_sav_biomes
#
# plot_biotic_df <- as.data.frame(plot_biotic_map, xy=T) %>%
#   filter(!is.na(layer)) %>%
#   mutate(layer = ifelse(layer == 0, NA, layer))
#
#
# p1 <- ggplot() +
#   geom_raster(data = plot_biotic_df, aes(fill = layer,
#                                          x=x, y=y)) +
#   scale_fill_viridis_c(limits = c(0, 1),
#                        values = c(0,0.3, 0.5, 0.6, 0.7, 1),
#                        direction = -1) +
#   coord_quickmap() +
#   theme_void() +
#   ylim(-60, 90) +
#   guides(fill = guide_colorbar(title="Biotic dispersal\nprevalence",
#                                barwidth = 12,
#                                title.vjust = 1.5,
#                                label.vjust = 1.6,
#                                barheight = 0.46)) +
#   theme(legend.position="bottom",
#         legend.title=element_text(size=10))
#
# my_inset <- function(){
#   par(mar = c(0,0,0,0),
#       mgp = c(1,0,0))
#   plot(NA,
#        xlim = c(-60, 80),
#        ylim = c(0, 1),
#        frame = F,
#        las = 1,
#        xlab = "",
#        ylab = "",
#        xaxt = "n",
#        yaxt = "n",
#        cex = 0.5)
#   axis(1, at = c(-60,0,60), labels = c("60°S", "0", "60°N"), las = 1, tck = -0.005,
#        cex.axis = 0.6)
#   axis(2, at = c(0,0.5,1), labels = c("0  ", "0.5  ", "1  "), las = 1, tck = -0.005,
#        cex.axis = 0.6)
#
#
#   set.seed(4)
#   points(biotic ~ jitter(latitude, 40),
#          data = both_mean,
#          pch = 16,
#          cex = 0.05,
#          col = "grey60")#rgb(0,0,0,0.1))
#
#   lines(plogis(gam_preds$latitude$predicted) ~ gam_preds$latitude$x,
#         col = "darkgreen")
#
#   polygon(y = plogis(c(gam_preds$latitude$conf.low,
#                        rev(gam_preds$latitude$conf.high))),
#           x = c(gam_preds$latitude$x,
#                 rev(gam_preds$latitude$x)),
#           col = rgb(34, 139, 34, maxColorValue = 255, alpha = 150), #2, 48, 32
#           border = NA)
# }
#
# ipak("grid")
# ipak("ggplotify")
#
# p2 <- as.grob(~my_inset())
#
# #grid.newpage()
# grid.draw(p1)
# vp = viewport(x=.15, y=.4, width=.15, height=.25)
# pushViewport(vp)
# grid.draw(p2)
# upViewport()
#
#
# dev.off()
