library('tidyverse')
library('stars')

# Set up working directories
data_wd <- '/nobackup1/efricke/VertNBS_data'
outputs_wd <- '/nobackup1/efricke/VertNBS_data'
raster_wd <- '/nobackup1/efricke/VertNBS_data'
raster_gee_wd <- '/nobackup1/efricke/VertNBS_data'
raster_outputs_wd <- '/nobackup1/efricke/VertNBS_visuals'
raster_intermediary_wd <- '/nobackup1/efricke/VertNBS_data/intermediary_rasters'

setwd(raster_outputs_wd)

# Get all the files for the current scenario
current_lo_files <- list.files(pattern = "c30_2020_current_trop_for_and_sav_biomes_lo_chunk")
current_hi_files <- list.files(pattern = "c30_2020_current_trop_for_and_sav_biomes_hi_chunk")

# Get all the files for diff
diff_lo_files <- list.files(pattern = "c30_2020_diff_trop_for_and_sav_biomes_lo_chunk")
diff_hi_files <- list.files(pattern = "c30_2020_diff_trop_for_and_sav_biomes_hi_chunk")


# Loop over these to build up the whole image.
# 1. current lo
c30_2020_current_trop_for_and_sav_biomes_lo <- read_stars(current_lo_files[1], proxy = F)
c30_2020_current_trop_for_and_sav_biomes_lo[[1]][is.na(c30_2020_current_trop_for_and_sav_biomes_lo[[1]])] <- 0

for(i in 2:length(current_lo_files)){
  chunk <- read_stars(current_lo_files[i], proxy = F)
  
  chunk[[1]][is.na(chunk[[1]])] <- 0
  c30_2020_current_trop_for_and_sav_biomes_lo[[1]] <- c30_2020_current_trop_for_and_sav_biomes_lo[[1]] + chunk[[1]]
  
  print(paste("current lo", i))
}

write_stars(c30_2020_current_trop_for_and_sav_biomes_lo,
            dsn = "c30_2020_current_trop_for_and_sav_biomes_lo_fullmap.tif")

rm(chunk, c30_2020_current_trop_for_and_sav_biomes_lo)
gc()

# 2. current hi
c30_2020_current_trop_for_and_sav_biomes_hi <- read_stars(current_hi_files[1], proxy = F)
c30_2020_current_trop_for_and_sav_biomes_hi[[1]][is.na(c30_2020_current_trop_for_and_sav_biomes_hi[[1]])] <- 0

for(i in 2:length(current_hi_files)){
  chunk <- read_stars(current_hi_files[i], proxy = F)
  
  chunk[[1]][is.na(chunk[[1]])] <- 0
  c30_2020_current_trop_for_and_sav_biomes_hi[[1]] <- c30_2020_current_trop_for_and_sav_biomes_hi[[1]] + chunk[[1]]
  
  print(paste("current hi", i))
}

write_stars(c30_2020_current_trop_for_and_sav_biomes_hi,
            dsn = "c30_2020_current_trop_for_and_sav_biomes_hi_fullmap.tif")

rm(chunk, c30_2020_current_trop_for_and_sav_biomes_hi)
gc()

# 3. diff lo
c30_2020_diff_trop_for_and_sav_biomes_lo <- read_stars(diff_lo_files[1], proxy = F)
c30_2020_diff_trop_for_and_sav_biomes_lo[[1]][is.na(c30_2020_diff_trop_for_and_sav_biomes_lo[[1]])] <- 0

for(i in 2:length(diff_lo_files)){
  chunk <- read_stars(diff_lo_files[i], proxy = F)
  
  chunk[[1]][is.na(chunk[[1]])] <- 0
  c30_2020_diff_trop_for_and_sav_biomes_lo[[1]] <- c30_2020_diff_trop_for_and_sav_biomes_lo[[1]] + chunk[[1]]
  
  print(paste("diff lo", i))
}

write_stars(c30_2020_diff_trop_for_and_sav_biomes_lo,
            dsn = "c30_2020_diff_trop_for_and_sav_biomes_lo_fullmap.tif")

rm(chunk, c30_2020_diff_trop_for_and_sav_biomes_lo)
gc()

# 4. diff hi
c30_2020_diff_trop_for_and_sav_biomes_hi <- read_stars(diff_hi_files[1], proxy = F)
c30_2020_diff_trop_for_and_sav_biomes_hi[[1]][is.na(c30_2020_diff_trop_for_and_sav_biomes_hi[[1]])] <- 0

for(i in 2:length(diff_hi_files)){
  chunk <- read_stars(diff_hi_files[i], proxy = F)
  
  chunk[[1]][is.na(chunk[[1]])] <- 0
  c30_2020_diff_trop_for_and_sav_biomes_hi[[1]] <- c30_2020_diff_trop_for_and_sav_biomes_hi[[1]] + chunk[[1]]
  
  print(paste("diff hi", i))
}

write_stars(c30_2020_diff_trop_for_and_sav_biomes_hi,
            dsn = "c30_2020_diff_trop_for_and_sav_biomes_hi_fullmap.tif")

rm(chunk, c30_2020_diff_trop_for_and_sav_biomes_hi)
gc()
