# Load packages ----------------------------------------------------------------
sapply(list.files("./R", full.names = T), source)

ipak(c("tidyverse"))


# Load move dataframe
move_df <- read.csv("./data/land_use_dependent_dispersal_data/tidy/move_df_spat.csv")


# Want to reconsider the way we handle the animal groupings
marine_foraging_birds <- c(
  "Alca torda",
  "Fratercula arctica",
  "Uria aalge",
  "Uria lomvia",
  "Hydroprogne caspia",
  "Ichthyaetus audouinii",
  "Larus marinus",
  "Onychoprion fuscatus",
  "Rissa tridactyla",
  "Sterna dougallii",
  "Pelecanus occidentalis",
  "Phaethon lepturus",
  "Phoebastria irrorata",
  "Thalassarche carteri",
  "Thalassarche chrysostoma",
  "Calonectris diomedea",
  "Puffinus puffinus",
  "Puffinus yelkouan",
  "Thalassoica antarctica",
  "Eudyptes robustus",
  "Pygoscelis adeliae",
  "Fregata aquila",
  "Fregata magnificens",
  "Morus bassanus",
  "Sula dactylatra")

marine_foraging_mammals <- c(
  "Halichoerus grypus",
  "Balaenoptera musculus",
  "Balaenoptera physalus")

raptor_orders <- c("Accipitriformes",
                   "Falconiformes",
                   "Strigiformes")

# Water and wading birds?

move_df %>%
  filter(class == "Aves",
         !individual.taxon.canonical.name %in% marine_foraging_birds,
         !order %in% raptor_orders) %>%
  dplyr::select(order, family, individual.taxon.canonical.name) %>%
  dplyr::arrange(order, family) %>%
  unique()

ww_bird_orders <- c("Anseriformes",
                    "Charadriiformes",
                    "Ciconiiformes",
                    "Gruiformes", # Be careful here. Currently just a crane and a coot, but there are land rails too...
                    "Pelecaniformes")
ww_bird_species <- c("Phalacrocorax carbo")

move_df$animal_group[which(move_df$individual.taxon.canonical.name %in% marine_foraging_birds)] <- "bird_marine"
# Raptors
move_df$animal_group[which(move_df$order %in% raptor_orders)] <- "bird_raptor"

# Marine mammals
move_df$animal_group[which(move_df$individual.taxon.canonical.name %in% marine_foraging_mammals)] <- "mamm_marine"

move_df$animal_group %>% table()



# Manipulate human footprint data so that the ocean and land
# footprint indices are at the same scale

# Note that the range of the ocean values is: 0, 11.7999 for the entire raster
# and the land max is ~50.

center_fun <- function(x){ x - mean(x, na.rm - T)}
range01 <- function(x){ (x - min(x, na.rm - T)) / (max(x, na.rm - T) - min(x, na.rm - T))}

# # Uncomment if doing both marine and terrestrial species (and NDVI and chlor)
# move_df <- move_df %>%
#   mutate(land_fp_n = land_fp / max(land_fp, na.rm = T),
#          ocean_fp_n = ocean_fp / max(ocean_fp, na.rm = T),
#          ndvi_terra_n = ndvi_terra / max(ndvi_terra, na.rm = T),
#          chlor_terra_n = (log(chlor_terra) + 2) / max((log(chlor_terra) + 2), na.rm = T)) %>%
#   mutate(fp = ifelse(!is.na(land_fp_n), land_fp_n, ocean_fp_n),
#          product = ifelse(!is.na(ndvi_terra_n), ndvi_terra_n, chlor_terra_n))

move_df <- move_df %>%
  mutate(land_fp_n = land_fp / max(land_fp, na.rm = T),
         ndvi_terra_n = ndvi_terra / max(ndvi_terra, na.rm = T)) %>%
  mutate(fp = land_fp_n,
         product = ndvi_terra_n)


# Get a mean_fp value and production value as well and join this to move_df
mean_fp <- move_df %>%
  group_by(frug_event_id) %>%
  dplyr::summarise(mean_fp = mean(fp, na.rm = T)) %>%
  dplyr::select(frug_event_id, mean_fp)

move_df <- move_df %>%
  left_join(mean_fp)

# Get a mean_fp value and production value as well and join this to move_df
mean_product <- move_df %>%
  group_by(frug_event_id) %>%
  dplyr::summarise(mean_product = mean(product, na.rm = T)) %>%
  dplyr::select(frug_event_id, mean_product)

move_df <- move_df %>%
  left_join(mean_product)


# We only want terrestrial data in here
terr_df <- move_df %>%
  filter(!animal_group %in% c("mamm_cetacean", "bird_marine"),
         individual.taxon.canonical.name != "Halichoerus grypus")

# Want to add some metadata
move_metadata <- read.csv("./data/land_use_dependent_dispersal_data/movebank_metadata.csv", header = T)

# Filter to ones that are included in terr_df
move_metadata <- move_metadata %>%
  filter(name %in% terr_df$study.id |
           id %in% terr_df$study.id) %>%
  mutate(study.id = ifelse(name %in% terr_df$study.id,
                           name, id))

# And it looks like there's some duplication where a study was added both
# by its "id" and by its "name" - will remove these
terr_df <- terr_df %>%
  filter(!study.id %in% c("193545363", "3109235", "492444603", "9480191"))

# Only want to keep ones with the CC_0 license
terr_df2 <- terr_df %>%
  left_join(move_metadata %>%
              dplyr::select(study.id, name, id,
                            license_type, license_terms)) %>%
  filter(license_type %in% c("CC_0"))

# And here is the model with random slopes and intercepts by species and intercepts by frugivory event
ipak("lme4")
terr_df2 <- terr_df2[complete.cases(terr_df2 %>%
                          dplyr::select(time_diff_min,
                                        mean_fp,
                                        mean_product)),]
terr_mod_mean <- lmer(log(displacement) ~ log(time_diff_min) + mean_fp + mean_product +
                        (log(time_diff_min)|individual.taxon.canonical.name) +
                        (1|frug_event_id),
                 REML = F,
                 control = lmerControl(optimizer = "Nelder_Mead"),
                 data = terr_df2)

summary(terr_mod_mean)

terr_mod_0 <- lmer(log(displacement) ~ log(time_diff_min) + mean_product +
                        (log(time_diff_min)|individual.taxon.canonical.name) +
                        (1|frug_event_id),
                      REML = F,
                      control = lmerControl(optimizer = "Nelder_Mead"),
                      data = terr_df2)
# Use a likelihood ratio test to assess statistical significance
anova(terr_mod_mean, terr_mod_0)

