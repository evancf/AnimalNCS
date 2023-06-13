source(list.files("./R", full.names = T))

ipak(c("tidyverse"))

# Load site covariate data -----------------------------------------------------


# Load coord_dat
coord_covars <- read.csv(file = "./data/coord_covars.csv", row.names = 1)


# Load mammal and bird presence by site data -----------------------------------

# Note that this code depends on an updated run of the code in
# "get_species_presence_by_sitecoords.R", which takes hours to run and
# depends on you having IUCN Red List / BirdLife International range
# maps downloaded.

# Load these in if possible (otherwise have to run the above mentioned file)
load("./data/spatial data/bird_presence_by_sitecoords_wide.RData")
load("./data/spatial data/mamm_presence_by_sitecoords_wide.RData")

# Make sure the column numbers match (this is n+1 number of sites)
dim(mamm_presence_by_sitecoords_wide)
dim(bird_presence_by_sitecoords_wide)

colnames(mamm_presence_by_sitecoords_wide) %>% head()
colnames(bird_presence_by_sitecoords_wide) %>% head()


mamm_presence_by_sitecoords_wide[1:6, 1:6]
bird_presence_by_sitecoords_wide[1:6, 1:6]
bird_presence_by_sitecoords_wide$SCINAME %>% unique() %>% length()

# Want to use "binomial" for scientific name rather than SCINAME
bird_presence_by_sitecoords_wide <- bird_presence_by_sitecoords_wide %>%
  rename(binomial = SCINAME)

(colnames(mamm_presence_by_sitecoords_wide) == colnames(bird_presence_by_sitecoords_wide)) %>% all()


# Load in datasets from Fricke et al. 2022 -------------------------------------
load(url("https://github.com/evancf/global-dispersal-change/blob/main/data-and-model-outputs/range.matrices.RData?raw=true"))
load(url("https://github.com/evancf/global-dispersal-change/blob/main/data-and-model-outputs/all.bird.traits.RData?raw=true"))
load(url("https://github.com/evancf/global-dispersal-change/blob/main/data-and-model-outputs/all.mamm.traits.RData?raw=true"))

# Also load predicted values for disperser species
load("./data/disperser_preds.RData")

# This is a big compilation of mammal traits, which we will also use
# https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3344
comb_dat <- read.csv("./data/COMBINE_archives/trait_data_imputed.csv")

# Need to deal with mammal names in IUCN being unequal to the ones in
# Phylacine and the Fricke et al. 2022 analysis.

# First, want to ignore subspecies
mamm_presence_by_sitecoords_wide <- mamm_presence_by_sitecoords_wide %>%
  mutate(binomial = word(binomial, 1, 2))

# Next want to figure out how the species are called in phylacine
phylacine_name <- mamm_presence_by_sitecoords_wide %>%
  dplyr::select(binomial) %>%
  left_join(comb_dat %>% dplyr::select(iucn2020_binomial,
                                       phylacine_binomial,
                                       order,
                                       family,
                                       genus,
                                       species),
            by = c("binomial" = "iucn2020_binomial"))

phylacine_name %>%
  filter(phylacine_binomial %>% is.na()) %>%
  pull("binomial")

gen <- "Baeodon"
sp <- "alleni"

bin <- "Cacajao hosom"
bin %in% spp_names_mamm

tibble(spp_names_mamm = spp_names_mamm) %>% filter(word(spp_names_mamm, 1) == gen)
tibble(spp_names_mamm = spp_names_mamm) %>% filter(word(spp_names_mamm, 2) == sp)

# Get the phylacine names
# First some that are just NAs
phyla_bin_na_inds <- which(is.na(phylacine_name$phylacine_binomial))
phylacine_name$phylacine_binomial[phyla_bin_na_inds] <- phylacine_name$binomial[phyla_bin_na_inds] %>%
  # These ones are just fine as is
  plyr::revalue(c("Capricornis milneedwardsii" = "Capricornis milneedwardsii",
                  "Capricornis thar" = "Capricornis thar",
                  "Fukomys mechowii" = NA, # No dispersal estimates anyway for this rodent anyway, so wont matter
                  "Hipposideros commersoni" = "Hipposideros commersoni",
                  "Hipposideros gigas" = "Hipposideros gigas",
                  "Hipposideros vittatus" = "Hipposideros vittatus",
                  "Lagothrix cana" = "Lagothrix cana",
                  "Lagothrix lugens" = "Lagothrix lugens",
                  "Megaderma lyra" = "Megaderma lyra",
                  "Naemorhedus griseus" = "Naemorhedus griseus",
                  "Ovis orientalis" = "Ovis orientalis",
                  "Pipistrellus ariel" = "Pipistrellus ariel",
                  "Pipistrellus musciculus" = "Pipistrellus musciculus",
                  "Pipistrellus pulveratus" = "Pipistrellus pulveratus",
                  "Saguinus fuscicollis" = "Saguinus fuscicollis",
                  "Sapajus macrocephalus" = "Sapajus macrocephalus",
                  "Tarsius bancanus" = "Tarsius bancanus"))

# Next ones that are listed as not recognized - basically a bunch of bats
phyla_bin_nr_inds <- which(phylacine_name$phylacine_binomial == "Not recognised")
phylacine_name$phylacine_binomial[phyla_bin_nr_inds] <- phylacine_name$binomial[phyla_bin_nr_inds] %>%
  # These ones are just fine as is
  plyr::revalue(c("Abrothrix manni" = NA, # Rodent
                  "Cacajao melanocephalus"  = "Cacajao hosomi",
                  "Canis lupaster"  = "Canis aureus",
                  "Cercopithecus pogonias"  = "Cercopithecus mona",
                  "Cercopithecus wolfi"  = "Cercopithecus mona",
                  "Drymoreomys albimaculatus"  = NA, # Rodent
                  "Hipposideros khasiana"  = "Hipposideros larvatus",
                  "Lonchophylla cadenai"  = "Lonchophylla thomasi", # Functionally closest, so will use these estimates. Recently described species
                  "Murina bicolor"  = "Murina leucogaster", # Similar
                  "Murina gracilis"  = "Murina leucogaster", # Similar
                  "Murina recondita"  = "Murina leucogaster", # Similar
                  "Myotis annatessae"  = "Myotis muricola",
                  "Myotis badius"  = "Myotis siligorensis",
                  "Myotis indochinensis"  = "Myotis montivagus",
                  "Myotis secundus"  = "Myotis yanbarensis",
                  "Plecotus ariel"  = "Plecotus austriacus",
                  "Plecotus homochrous"  = "Plecotus auritus",
                  "Rhinolophus xinanzhongguoensis" = "Rhinolophus ferrumequinum",
                  "Rhipidomys itoan"  = NA, # Rodent
                  "Scotophilus trujilloi"  = "Scotophilus dinganii",
                  "Sturnira bakeri" = "Sturnira lilium"))


filter(phylacine_name, phylacine_binomial == "Not recognised") %>%
  pull(binomial)

# Okay, clean this up a little
phylacine_name <- phylacine_name[!duplicated(phylacine_name$binomial),]

# And switch out the IUCN names for the phylacine names
mamm_presence_by_sitecoords_wide <- mamm_presence_by_sitecoords_wide %>%
  left_join(phylacine_name %>% dplyr::select(binomial, phylacine_binomial)) %>%
  mutate(binomial = phylacine_binomial) %>%
  dplyr::select(-phylacine_binomial)




# Want to figure out which birds and mammals should be removed from analysis
# For example, humans, marine organisms

# Determine which species to include in analyses
# Skip marine mammals
marine.mamm.to.skip <- filter(all.mamm.traits, marine == 1)$name.m

not.terrestrial <- matrix(colSums(m.mamm.current[-which(rownames(m.mamm.current) %in% marine.mamm.to.skip),]),
                          nrow = 142,
                          ncol = 360, byrow = T) == 0

# Only focus on mammal orders including frugivores (note that this
# overlooks a role for rodents and other ecologically similar small mammals
# in dispersal of fleshy-fruited plants)
orders.bats <- "Chiroptera"
orders.primates <- "Primates"
orders.mamm.carns <- c("Carnivora", "Cingulata", "Dasyuromorphia",
                       "Didelphimorphia", "Eulipotyphla", "Microbiotheria",
                       "Monotremata", "Paucituberculata", "Pholidota",
                       "Tubulidentata")
orders.mamm.herbs <- c("Artiodactyla", "Cetartiodactyla", "Dermoptera",  # Note that the traits dataframe and net.long use different order names...
                       "Diprotodontia", "Hyracoidea", "Litopterna",
                       "Notoungulata", "Perissodactyla", "Pilosa",
                       "Proboscidea")

orders.mamm.to.skip <- filter(all.mamm.traits, !Order.1.2 %in% c(orders.bats, orders.primates,
                                                                 orders.mamm.carns, orders.mamm.herbs) )$name.m

spp.mamm.to.skip <- c(marine.mamm.to.skip, orders.mamm.to.skip, "Homo sapiens")

# Remove these species from consideration
mamm_presence_by_sitecoords_wide <- mamm_presence_by_sitecoords_wide[!(mamm_presence_by_sitecoords_wide$binomial %in% spp.mamm.to.skip), ]

# There are a few NAs as binomial - these are just things that need to be omitted
mamm_presence_by_sitecoords_wide$binomial %>% is.na() %>% table()
mamm_presence_by_sitecoords_wide <- mamm_presence_by_sitecoords_wide[!is.na(mamm_presence_by_sitecoords_wide$binomial),]


# Next need to reconcile synonyms so they aren't double counted
mamm_presence_by_sitecoords_wide$binomial %>% length()
mamm_presence_by_sitecoords_wide$binomial %>% unique() %>% length()

mamm_presence_by_sitecoords_wide <- mamm_presence_by_sitecoords_wide %>%
  group_by(binomial) %>%
  summarise_each(sum) %>%
  mutate_if(is.numeric, ~ as.logical(1 * (. > 0))) %>%
  ungroup()



# Prepare for various dispersal estimates --------------------------------------

# Develop dispersal estimates based on a modification of the
# effective seed dispersal modeling in Fricke et al. 2022 Science

source("./analysis/disp_estimates_Fricke2022.R")

# Assign distance threshold
dist_thresh_mamm <- 1000
dist_thresh_bird <- 1000

# Effective dispersal
# Will first add ones with landscape effect
disp_ests_eff_le_mean <- dispersal_estimates(dist_threshold_mamm = dist_thresh_mamm,
                            dist_threshold_bird = dist_thresh_bird,
                            include_landscape_effect = T,
                            sum_or_mean = "mean",
                            type = "effective_dispersal")

disp_covars <- coord_covars %>%
  dplyr::select(coord_id) %>%
  mutate(disp_mamm_eff_le_mean = disp_ests_eff_le_mean[["mamm"]] %>% unlist(),
         disp_bird_eff_le_mean = disp_ests_eff_le_mean[["bird"]] %>% unlist(),
         disp_both_eff_le_mean = disp_ests_eff_le_mean[["both"]] %>% unlist())

# Will then add ones without landscape effect
disp_ests_eff_no_le_mean <- dispersal_estimates(dist_threshold_mamm = dist_thresh_mamm,
                                             dist_threshold_bird = dist_thresh_bird,
                                             include_landscape_effect = F,
                                             sum_or_mean = "mean",
                                             type = "effective_dispersal")

disp_covars <- disp_covars %>%
  mutate(disp_mamm_eff_no_le_mean = disp_ests_eff_no_le_mean[["mamm"]] %>% unlist(),
         disp_bird_eff_no_le_mean = disp_ests_eff_no_le_mean[["bird"]] %>% unlist(),
         disp_both_eff_no_le_mean = disp_ests_eff_no_le_mean[["both"]] %>% unlist())

# Participation
# Will first add ones with landscape effect
disp_ests_part_le_mean <- dispersal_estimates(dist_threshold_mamm = dist_thresh_mamm,
                                             dist_threshold_bird = dist_thresh_bird,
                                             include_landscape_effect = T,
                                             sum_or_mean = "mean",
                                             type = "participation")

disp_covars <- disp_covars %>%
  mutate(disp_mamm_part_le_mean = disp_ests_part_le_mean[["mamm"]] %>% unlist(),
         disp_bird_part_le_mean = disp_ests_part_le_mean[["bird"]] %>% unlist(),
         disp_both_part_le_mean = disp_ests_part_le_mean[["both"]] %>% unlist())

# Will then add ones without landscape effect
disp_ests_part_no_le_mean <- dispersal_estimates(dist_threshold_mamm = dist_thresh_mamm,
                                                dist_threshold_bird = dist_thresh_bird,
                                                include_landscape_effect = F,
                                                sum_or_mean = "mean",
                                                type = "participation")

disp_covars <- disp_covars %>%
  mutate(disp_mamm_part_no_le_mean = disp_ests_part_no_le_mean[["mamm"]] %>% unlist(),
         disp_bird_part_no_le_mean = disp_ests_part_no_le_mean[["bird"]] %>% unlist(),
         disp_both_part_no_le_mean = disp_ests_part_no_le_mean[["both"]] %>% unlist())

# The sum version
# Effective dispersal
# Will first add ones with landscape effect
disp_ests_eff_le_sum <- dispersal_estimates(dist_threshold_mamm = dist_thresh_mamm,
                                             dist_threshold_bird = dist_thresh_bird,
                                             include_landscape_effect = T,
                                             sum_or_mean = "sum",
                                             type = "effective_dispersal")

disp_covars <- disp_covars %>%
  mutate(disp_mamm_eff_le_sum = disp_ests_eff_le_sum[["mamm"]] %>% unlist(),
         disp_bird_eff_le_sum = disp_ests_eff_le_sum[["bird"]] %>% unlist(),
         disp_both_eff_le_sum = disp_ests_eff_le_sum[["both"]] %>% unlist())

# Will then add ones without landscape effect
disp_ests_eff_no_le_sum <- dispersal_estimates(dist_threshold_mamm = dist_thresh_mamm,
                                                dist_threshold_bird = dist_thresh_bird,
                                                include_landscape_effect = F,
                                                sum_or_mean = "sum",
                                                type = "effective_dispersal")

disp_covars <- disp_covars %>%
  mutate(disp_mamm_eff_no_le_sum = disp_ests_eff_no_le_sum[["mamm"]] %>% unlist(),
         disp_bird_eff_no_le_sum = disp_ests_eff_no_le_sum[["bird"]] %>% unlist(),
         disp_both_eff_no_le_sum = disp_ests_eff_no_le_sum[["both"]] %>% unlist())

# Participation
# Will first add ones with landscape effect
disp_ests_part_le_sum <- dispersal_estimates(dist_threshold_mamm = dist_thresh_mamm,
                                              dist_threshold_bird = dist_thresh_bird,
                                              include_landscape_effect = T,
                                              sum_or_mean = "sum",
                                              type = "participation")

disp_covars <- disp_covars %>%
  mutate(disp_mamm_part_le_sum = disp_ests_part_le_sum[["mamm"]] %>% unlist(),
         disp_bird_part_le_sum = disp_ests_part_le_sum[["bird"]] %>% unlist(),
         disp_both_part_le_sum = disp_ests_part_le_sum[["both"]] %>% unlist())

# Will then add ones without landscape effect
disp_ests_part_no_le_sum <- dispersal_estimates(dist_threshold_mamm = dist_thresh_mamm,
                                                 dist_threshold_bird = dist_thresh_bird,
                                                 include_landscape_effect = F,
                                                 sum_or_mean = "sum",
                                                 type = "participation")

disp_covars <- disp_covars %>%
  mutate(disp_mamm_part_no_le_sum = disp_ests_part_no_le_sum[["mamm"]] %>% unlist(),
         disp_bird_part_no_le_sum = disp_ests_part_no_le_sum[["bird"]] %>% unlist(),
         disp_both_part_no_le_sum = disp_ests_part_no_le_sum[["both"]] %>% unlist())



# Write the dispersal part as as csv -------------------------------------------
write.csv(file = "./data/disp_covars.csv", disp_covars)

