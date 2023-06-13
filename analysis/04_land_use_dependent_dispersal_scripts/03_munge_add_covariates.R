# Load packages ----------------------------------------------------------------
sapply(list.files("./R", full.names = T), source)

ipak(c("tidyverse",
       "move",
       "lutz",
       "lubridate"))


# Read in data
manipulated_files <- list.files("./data/land_use_dependent_dispersal_data/tidy/Movebank manipulated", full.names = T)
move_list <- list()
for(i in manipulated_files){
  move_list[[i]] <- read.csv(i, header = T, row.names = 1)
  if("study.name" %in% colnames(move_list[[i]])){
    colnames(move_list[[i]])[which(colnames(move_list[[i]]) == "study.name")] <- "study.id"

    colnames(move_list[[i]])[which(colnames(move_list[[i]]) == "sensor.type")] <- "sensor.type.id"
  }
}

move_df <- do.call(rbind, move_list)

dim(move_df)
table(move_df$individual.taxon.canonical.name) %>% sort()
table(move_df$frug_event_id) %>% hist()

# Fix up some species name issues
move_df[which(move_df$study.id == "2760899"),
        "individual.taxon.canonical.name"] <- "Fregata magnificens"
move_df[which(move_df$study.id == "1326534946"),
        "individual.taxon.canonical.name"] <- "Rissa tridactyla"
move_df[which(move_df$study.id == "113636867"),
        "individual.taxon.canonical.name"] <- "Rissa tridactyla"



# Decide which taxa to focus on
move_spp <- unique(move_df$individual.taxon.canonical.name)
bad_spp <- move_spp[!move_spp %in% c(traitdata::elton_birds$scientificNameStd,
                         traitdata::elton_mammals$scientificNameStd)]

move_df %>% filter(individual.taxon.canonical.name %in% bad_spp) %>%
  filter(!duplicated(study.id))


# Adding in covariates ---------------------------------------------------------

# Traits

etb <- traitdata::elton_birds
etm <- traitdata::elton_mammals

# Fix up a couple things with these..
library(taxizedb)

if(!"etm_orders.csv" %in% list.files("./data/land_use_dependent_dispersal_data/")){
  etm_orders <- taxize::tax_name(unique(etm$Family), get = "Order", db = "itis")
  etm_orders$order[which(etm_orders$query == "Ziphiidae")] <- "Cetacea"
  write.csv(etm_orders, file = "./data/land_use_dependent_dispersal_data/etm_orders.csv")
} else{
  etm_orders <- read.csv("./data/land_use_dependent_dispersal_data/etm_orders.csv", header = T, row.names = 1)
}

# Fix this one up
etb$scientificNameStd[which(etb$scientificNameStd == "Sterna fuscata")] <- "Onychoprion fuscatus"


etb <- etb %>%
  dplyr::select("scientificNameStd", "Order", "Family", "Genus",
                "BodyMass.Value") %>%
  rename("individual.taxon.canonical.name" = "scientificNameStd",
         "order" = "Order",
         "family" = "Family",
         "genus" = "Genus",
         "mass" = "BodyMass.Value") %>%
  mutate(class = "Aves")

etm <- etm %>%
  left_join(etm_orders %>% dplyr::select(query, order), by = c("Family" = "query")) %>%
  dplyr::select("scientificNameStd", "order", "Family", "Genus",
                "BodyMass.Value")  %>%
  rename("individual.taxon.canonical.name" = "scientificNameStd",
         "family" = "Family",
         "genus" = "Genus",
         "mass" = "BodyMass.Value") %>%
  mutate(class = "Mammalia")

et <- bind_rows(etb, etm)

trait_summary <- tibble(individual.taxon.canonical.name = move_df$individual.taxon.canonical.name %>% unique() %>% sort()) %>%
  left_join(et)



# Want to include an animal group term
orders_bats <- "Chiroptera"
orders_primates <- "Primates"
orders_mamm_carns <- c("Carnivora", "Cingulata", "Dasyuromorphia",
                       "Didelphimorphia", "Eulipotyphla", "Microbiotheria",
                       "Monotremata", "Paucituberculata", "Pholidota",
                       "Tubulidentata")
orders_mamm_herbs <- c("Artiodactyla", "Cetartiodactyla", "Dermoptera",  # Note that the traits dataframe and net_long use different order names___
                       "Diprotodontia", "Hyracoidea", "Litopterna",
                       "Notoungulata", "Perissodactyla", "Pilosa",
                       "Proboscidea", "Lagomorpha")

orders_mamm_cetacean <- "Cetacea"

trait_summary <- trait_summary %>%
  mutate(animal_group = ifelse(class == "Aves", "birds", NA)) %>%
  mutate(animal_group = ifelse(order %in% orders_bats, "bats", animal_group)) %>%
  mutate(animal_group = ifelse(order %in% orders_primates, "primates", animal_group)) %>%
  mutate(animal_group = ifelse(order %in% orders_mamm_carns, "mamm_carns", animal_group)) %>%
  mutate(animal_group = ifelse(order %in% orders_mamm_herbs, "mamm_herbs", animal_group)) %>%
  mutate(animal_group = ifelse(order %in% orders_mamm_cetacean, "mamm_cetacean", animal_group))


# Join the animal group info to the movement dataframe
move_df <- move_df %>% left_join(trait_summary)


# Some data manipulation -------------------------------------------------------

# We're going to remove hypothetical frugivory events where there's
# apparently no movement. These could be fallen tags or just animals
# that aren't active (so likely won't be eating fruits)

frug_event_ids_to_remove <- move_df %>%
  group_by(frug_event_id) %>%
  summarise(disp_mean = mean(displacement),
            disp_var = var(displacement)) %>%
  filter(disp_mean == 0 & disp_var == 0) %>%
  dplyr::select(frug_event_id) %>%
  pull(frug_event_id)

# Also for analysis, we don't need time = 0 and distance = 0
move_df <- move_df %>% filter(time_diff_min > 0) %>%
  filter(!frug_event_id %in% frug_event_ids_to_remove)

# In order to log-transform, we will change from displacement < 1 to 1
move_df$displacement[which(move_df$displacement < 1)] <- 1 #min(move_df$displacement[which(move_df$displacement != 0)])

# Remove distance outliers
move_df$displacement %>% log() %>% hist(breaks = 60)
abline(v = log(2000000))
move_df <- move_df %>% filter(displacement < 2000000)


# Will just use short term sensors here. Mainly GPS and radio transmitters
short_term_sensors <- c(653,673,2365682,1239574236,
                        "gps", "radio-transmitter", "natural-mark", "acoustic-tellemetry")
move_df <- move_df %>%
  filter(sensor.type.id %in% short_term_sensors)



# Want to remove a few species that have omove_df values
move_df <- move_df %>%
  filter(!individual.taxon.canonical.name %in% c("Calidris pusilla",
                                                 "Prunella modularis",
                                                 "Sylvia atricapilla",
                                                 "Sylvia borin"))

# Will also remove some species that are human-associated
move_df <- move_df %>%
  filter(!individual.taxon.canonical.name %in%
           c("Bos taurus",
             "Columba livia",
             "Felis catus",
             "Homo sapiens",
             "Ovis aries"))

# There are a few where there isn't species-level data, so we will
# just remove these
move_df <- move_df[complete.cases(move_df), ]



# Lastly, ensure that things are arranged correctly so that the
# factor levels make sense

move_df <- move_df %>%
  arrange(individual.taxon.canonical.name)

write.csv(move_df, file = "./data/tidy/move_df.csv")
# move_df <- read.csv("./data/tidy/move_df.csv")


