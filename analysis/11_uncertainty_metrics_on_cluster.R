# Choose the working directory to put all these scripts. Will just put these in
# a temp location because I'll upload them for execution on our slurm cluster

setwd("~/Downloads/")

for(i in 1:25){

  script_name <- paste0("uncertainty_metrics_on_cluster", i, ".R")

  sink(script_name)
  cat(paste0("

# Create an identifier for the filenames output by the script. This function
# will make a suffix particular to this output script


breaks <- matrix(c(1, 7,
                   8, 15,
                   16, 22,
                   23, 29,
                   30, 36,
                   37, 43,
                   44, 50,
                   51, 57,
                   58, 64,
                   65, 71,
                   72, 78,
                   79, 85,
                   86, 92,
                   93, 99,
                   100, 106,
                   107, 113,
                   114, 120,
                   121, 127,
                   128, 134,
                   135, 141,
                   142, 148,
                   149, 155,
                   156, 162,
                   163, 169,
                   170, 176),
                 ncol = 2,
                 byrow = T)

tile_indices <- breaks[", i ,", 1]:breaks[", i ,", 2]


library('tidyverse')
library('stars')

# Set up working directories
data_wd <- '/nobackup1/efricke/VertNBS_data'
outputs_wd <- '/nobackup1/efricke/VertNBS_data'
raster_wd <- '/nobackup1/efricke/VertNBS_data'
raster_gee_wd <- '/nobackup1/efricke/VertNBS_data'
raster_outputs_wd <- '/nobackup1/efricke/VertNBS_visuals'
raster_intermediary_wd <- '/nobackup1/efricke/VertNBS_data/intermediary_rasters'


suffix <- ", i ,"

create_file_suffix <- function(x, id = suffix){
  x <- gsub('.tif', '', x, fixed = T)
  x <- paste0(x, '_chunk', id, '.tif')
  x
}

#
setwd(raster_outputs_wd)

tif <- 'c30_2020_current_trop_for_and_sav_biomes.tif'
c30_2020_current_trop_for_and_sav_biomes <- read_stars(tif, proxy = TRUE)

current_files <- list.files(pattern = 'c30_2020_current_trop_for_and_sav_biomes_samp')
diff_files <- list.files(pattern = 'c30_2020_diff_trop_for_and_sav_biomes_samp')

tile_dim <- 1024
tiles <- st_tile(nrow(c30_2020_current_trop_for_and_sav_biomes),
                 ncol(c30_2020_current_trop_for_and_sav_biomes),
                 tile_dim, tile_dim)

# Check which tiles are all NA
na_tiles <- c()
for(i in seq_len(nrow(tiles))){
  tile <- read_stars(tif,
                     proxy = FALSE,
                     RasterIO = tiles[i, ])
  if(sum(!is.na(tile[[1]])) == 0){
    na_tiles <- c(na_tiles, i)
  }
  print(signif(i/nrow(tiles), 3))
}

# Remove the ones that are all NA
tiles <- tiles[-na_tiles,]

# Set up an output raster for the current uncertainty
c30_2020_current_trop_for_and_sav_biomes_lo <- read_stars(tif, proxy = F)
c30_2020_current_trop_for_and_sav_biomes_lo[[1]] <- NA
c30_2020_current_trop_for_and_sav_biomes_hi <- read_stars(tif, proxy = F)
c30_2020_current_trop_for_and_sav_biomes_hi[[1]] <- NA
gc()

for (i in seq_len(nrow(tiles))[tile_indices]) {
  tile <- read_stars(current_files,
                     proxy = FALSE,
                     RasterIO = tiles[i, ],
                     along = 'band')

  tile_lo <- st_apply(tile, 1:2, function(x) quantile(x, 0.025, na.rm = T))
  tile_hi <- st_apply(tile, 1:2, function(x) quantile(x, 0.975, na.rm = T))

  # Now place this onto the output file
  c30_2020_current_trop_for_and_sav_biomes_lo[[1]][tiles[i,1]:(tiles[i,1] + tiles[i,3] - 1),
                                                   tiles[i,2]:(tiles[i,2] + tiles[i,4] - 1)] <- tile_lo[[1]]

  # Now place this onto the output file
  c30_2020_current_trop_for_and_sav_biomes_hi[[1]][tiles[i,1]:(tiles[i,1] + tiles[i,3] - 1),
                                                   tiles[i,2]:(tiles[i,2] + tiles[i,4] - 1)] <- tile_hi[[1]]
}

write_stars(c30_2020_current_trop_for_and_sav_biomes_hi,
  dsn = 'c30_2020_current_trop_for_and_sav_biomes_hi.tif' %>% create_file_suffix())

write_stars(c30_2020_current_trop_for_and_sav_biomes_lo,
  dsn = 'c30_2020_current_trop_for_and_sav_biomes_lo.tif' %>% create_file_suffix())

rm(c30_2020_current_trop_for_and_sav_biomes_hi,
   c30_2020_current_trop_for_and_sav_biomes_lo)

gc()


# Set up an output raster for the diff uncertainty
tif <- 'c30_2020_diff_trop_for_and_sav_biomes.tif'
c30_2020_diff_trop_for_and_sav_biomes_lo <- read_stars(tif, proxy = F)
c30_2020_diff_trop_for_and_sav_biomes_lo[[1]] <- NA
c30_2020_diff_trop_for_and_sav_biomes_hi <- read_stars(tif, proxy = F)
c30_2020_diff_trop_for_and_sav_biomes_hi[[1]] <- NA
gc()

for (i in seq_len(nrow(tiles))[tile_indices]) {
  tile <- read_stars(diff_files,
                     proxy = FALSE,
                     RasterIO = tiles[i, ],
                     along = 'band')

  tile_lo <- st_apply(tile, 1:2, function(x) quantile(x, 0.025, na.rm = T))
  tile_hi <- st_apply(tile, 1:2, function(x) quantile(x, 0.975, na.rm = T))

  # Now place this onto the output file
  c30_2020_diff_trop_for_and_sav_biomes_lo[[1]][tiles[i,1]:(tiles[i,1] + tiles[i,3] - 1),
                                                   tiles[i,2]:(tiles[i,2] + tiles[i,4] - 1)] <- tile_lo[[1]]

  # Now place this onto the output file
  c30_2020_diff_trop_for_and_sav_biomes_hi[[1]][tiles[i,1]:(tiles[i,1] + tiles[i,3] - 1),
                                                   tiles[i,2]:(tiles[i,2] + tiles[i,4] - 1)] <- tile_hi[[1]]
}

write_stars(c30_2020_diff_trop_for_and_sav_biomes_hi,
  dsn = 'c30_2020_diff_trop_for_and_sav_biomes_hi.tif' %>% create_file_suffix())

write_stars(c30_2020_diff_trop_for_and_sav_biomes_lo,
  dsn = 'c30_2020_diff_trop_for_and_sav_biomes_lo.tif' %>% create_file_suffix())



"),fill=TRUE)
  sink()

  submit_script_name <- paste0("submit.uncertainty.metrics", i)
  sink(submit_script_name)
  cat(paste0("#!/bin/bash
#SBATCH -n 16 #request cpus
#SBATCH -N 1 #request nodes
#SBATCH -t 0-03:00 #time
#SBATCH -C centos7 #only centos7 nodes
#SBATCH -p sched_mit_hill #Run on sched_engaging_default partition
#SBATCH --mem=0 #
#SBATCH -o output_%j.txt #the job output will be directed to the file output_JOBID.txt
#SBATCH -e error_%j.txt #Errors to this
#SBATCH --mail-type=BEGIN,END #Mail at beginning and end of job
#SBATCH --mail-user=efricke@mit.edu
module load R/4.1.0
module add proj/7.2.1
module add gdal/3.2.1
module add geos/3.8.3
module add udunits2/2.2.26
module add engaging/cmake/3.5.2
module add gsl/2.5
module add gcc/9.3.0
R CMD BATCH /nobackup1/efricke/VertNBS_scripts/", script_name),
      fill=TRUE)
  sink()

}

