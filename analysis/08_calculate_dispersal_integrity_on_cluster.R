# The goal here is to calculate dispersal integrity for each pixel of
# a global map. Because this uses large rasters and large vector data on
# species ranges, it's not straightforward to perform these calcualtions
# in a memmory efficient way, so this breaks it up into chunks.

# This script that makes a bunch of individual scripts that can be run on
# a cluster because this is really memory/time intensive. This is a little
# clunky and I'm sure there's a more straightforward way to run this in parallel

# Choose the working directory to put all these scripts. Will just put these in
# a temp location because I'll upload them to our slurm cluster

setwd("~/Downloads/")

# Also, see below for the SBATCH submit scripts that are also produced.

# Further, we will duplicate this whole thing to get dispersal integrity
# estimates for 2020 below


for(i in 1:22){

  script_name <- paste0("calculate_dispersal_integrity_on_cluster_", i, ".R")

  sink(script_name)
  cat(paste0("

# This is the master script index used to perform this all in parallel (in a crude way)
script_index <- ", i ,"
n_scripts <- 22

# Here are the indices used later on
ind_seq <- seq(1, 41273, length.out = n_scripts + 1)
ind_start <- ind_seq[script_index]
ind_end <- ind_seq[script_index + 1]

# Output filename
output_filename <- paste0('out_list2000', script_index, '.rds')

# Directories
output_dir <- '/nobackup1/efricke/VertNBS_outputs'
data_dir <- '/nobackup1/efricke/VertNBS_data'
script_dir <- '/nobackup1/efricke/VertNBS_scripts'

time_1 <- Sys.time()
time_1

# Packages
library('tidyverse')
library('Rfast')
library('parallel')
library('sf')
library('stars')
library('tictoc')
library('raster')

# Some functions and objects ------------------

# Heres the equal area projection we will use
proj_crs <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'

# Prepare a few objects that will be used for making the dispersal predictions
pgerm <- function(yi, prob = 0.4){
  (exp(yi) * prob) / (1 + exp(yi) * prob - prob)
}

# Assign a threshold for long-distance dispersal.
# We will use 1000 meters
ldd_dist_mamm <- 1000
ldd_dist_bird <- 1000
my_sdlog <- 0.6 # Based on a meta analysis


# Load data ------------------
setwd(data_dir)

data_files <- list.files(full.names = T)

# Load fp and gfc rasters as stars objects but not as proxies
fp_stars <- read_stars('fp.tif',
                       proxy = F)
gfc_stars <- read_stars('gfc.tif',
                        proxy = F)

# Will use this raster to help make some conversions between rows/columns
# and cells in the stars object later on
fp_rasterlayer <- as(fp_stars, 'Raster')


# Load range data
mamm_range_filename <- data_files[grepl('mamm_ranges', data_files)]
bird_range_filename <- data_files[grepl('bird_ranges', data_files)]

mamm_ranges_sf <- st_read(mamm_range_filename)
bird_ranges_sf <- st_read(bird_range_filename)


# Load dispersal prediction dataframes
load('./preds_dfs.RData')

# And make sure the order is also exactly the same!
(preds_df_mamm$n_preds == mamm_ranges_sf$n_preds) %>% all() # If true, it is good to go

# And make sure the order is also exactly the same!
(preds_df_bird$n_preds == bird_ranges_sf$n_preds) %>% all() # If true, it is good to go




# Now going to make the raster of dispersal reduction ---------------------------
hr_x <- dim(gfc_stars)[1]
hr_y <- dim(gfc_stars)[2]

# Output the dispersal integrity values
set.seed(5)
ind_vals <- 1:hr_x

bbox_list <- mclapply(1:length(ind_vals), function(x) {
  st_bbox(extent(fp_rasterlayer,
                 1,
                 hr_y,
                 ind_vals[x],
                 ind_vals[x]),
          crs = proj_crs)
})


time_2 <- Sys.time()
time_2
out_list <- mclapply((ind_start:ind_end), function(x) {
  rc <- st_crop(fp_stars, bbox_list[[x]])
  rc_sf <- sf::st_as_sf(rc, na.rm = F, as_points = T)

  mamm_range_mat <- st_intersects(mamm_ranges_sf, rc_sf) %>% as.matrix()

  # get the land-use dependent probability of long distance dispersal at each cell
  p_ldd_mamm <- 1-plnorm(q = ldd_dist_mamm,
                         meanlog = log(preds_df_mamm$mamm_meanlog_inside %*%
                                         t(exp(-1.86 * fp_stars[[1]][ind_vals[x],]/50))),
                         sdlog = my_sdlog)

  # Get the product
  cell_prod_mamm <- colsums(t(plogis(-0.99 + 1.31 * gfc_stars[[1]][ind_vals[x],]/100) * t(p_ldd_mamm * mamm_range_mat)))
  cell_prod0_mamm <- colsums(preds_df_mamm$p_ldd0_mamm * mamm_range_mat)



  bird_range_mat <- st_intersects(bird_ranges_sf, rc_sf) %>% as.matrix()

  # get the land-use dependent probability of long distance dispersal at each cell
  p_ldd_bird <- 1-plnorm(q = ldd_dist_bird,
                         meanlog = log(preds_df_bird$bird_meanlog_inside %*%
                                         t(exp(-1.86 * fp_stars[[1]][ind_vals[x],]/50))),
                         sdlog = my_sdlog)

  # Get the product
  cell_prod_bird <- colsums(t(plogis(-0.99 + 1.31 * gfc_stars[[1]][ind_vals[x],]/100) * t(p_ldd_bird * bird_range_mat)))
  cell_prod0_bird <- colsums(preds_df_bird$p_ldd0_bird * bird_range_mat)

  (cell_prod_mamm + cell_prod_bird) / (cell_prod0_mamm + cell_prod0_bird)

})
time_3 <- Sys.time()
time_3
time_3 - time_2


# Save this file
setwd(output_dir)
saveRDS(out_list, file = output_filename)


"),fill=TRUE)
  sink()

  submit_script_name <- paste0("submit.script", i)
  sink(submit_script_name)
  cat(paste0("#!/bin/bash
#SBATCH -n 16 #request cpus
#SBATCH -N 1 #request nodes
#SBATCH -t 0-08:00 #time
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




# And now the 2020 version

for(i in 1:22){

  script_name <- paste0("calculate_dispersal_integrity2020_on_cluster_", i, ".R")

  sink(script_name)
  cat(paste0("

# This is the master script index used to perform this all in parallel (in a crude way)
script_index <- ", i ,"
n_scripts <- 22

# Here are the indices used later on
ind_seq <- seq(1, 41273, length.out = n_scripts + 1)
ind_start <- ind_seq[script_index]
ind_end <- ind_seq[script_index + 1]

# Output filename
output_filename <- paste0('out_list2020', script_index, '.rds')

# Directories
output_dir <- '/nobackup1/efricke/VertNBS_outputs'
data_dir <- '/nobackup1/efricke/VertNBS_data'
script_dir <- '/nobackup1/efricke/VertNBS_scripts'

time_1 <- Sys.time()
time_1

# Packages
library('tidyverse')
library('Rfast')
library('parallel')
library('sf')
library('stars')
library('tictoc')
library('raster')

# Some functions and objects ------------------

# Heres the equal area projection we will use
proj_crs <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'

# Prepare a few objects that will be used for making the dispersal predictions
pgerm <- function(yi, prob = 0.4){
  (exp(yi) * prob) / (1 + exp(yi) * prob - prob)
}

# Assign a threshold for long-distance dispersal.
# We will use 1000 meters
ldd_dist_mamm <- 1000
ldd_dist_bird <- 1000
my_sdlog <- 0.6 # Based on a meta analysis


# Load data ------------------
setwd(data_dir)

data_files <- list.files(full.names = T)

# Load fp and gfc rasters as stars objects but not as proxies
fp_stars <- read_stars('fp.tif',
                       proxy = F)
gfc_stars <- read_stars('gfc2020.tif',
                        proxy = F)

# Will use this raster to help make some conversions between rows/columns
# and cells in the stars object later on
fp_rasterlayer <- as(fp_stars, 'Raster')


# Load range data
mamm_range_filename <- data_files[grepl('mamm_ranges', data_files)]
bird_range_filename <- data_files[grepl('bird_ranges', data_files)]

mamm_ranges_sf <- st_read(mamm_range_filename)
bird_ranges_sf <- st_read(bird_range_filename)


# Load dispersal prediction dataframes
load('./preds_dfs.RData')

# And make sure the order is also exactly the same!
(preds_df_mamm$n_preds == mamm_ranges_sf$n_preds) %>% all() # If true, it is good to go

# And make sure the order is also exactly the same!
(preds_df_bird$n_preds == bird_ranges_sf$n_preds) %>% all() # If true, it is good to go




# Now going to make the raster of dispersal reduction ---------------------------
hr_x <- dim(gfc_stars)[1]
hr_y <- dim(gfc_stars)[2]

# Output the dispersal integrity values
set.seed(5)
ind_vals <- 1:hr_x

bbox_list <- mclapply(1:length(ind_vals), function(x) {
  st_bbox(extent(fp_rasterlayer,
                 1,
                 hr_y,
                 ind_vals[x],
                 ind_vals[x]),
          crs = proj_crs)
})


time_2 <- Sys.time()
time_2
out_list <- mclapply((ind_start:ind_end), function(x) {
  rc <- st_crop(fp_stars, bbox_list[[x]])
  rc_sf <- sf::st_as_sf(rc, na.rm = F, as_points = T)

  mamm_range_mat <- st_intersects(mamm_ranges_sf, rc_sf) %>% as.matrix()

  # get the land-use dependent probability of long distance dispersal at each cell
  p_ldd_mamm <- 1-plnorm(q = ldd_dist_mamm,
                         meanlog = log(preds_df_mamm$mamm_meanlog_inside %*%
                                         t(exp(-1.86 * fp_stars[[1]][ind_vals[x],]/50))),
                         sdlog = my_sdlog)

  # Get the product # This is a little clunky, but need to remember not to divide this raster by 100
  cell_prod_mamm <- colsums(t(plogis(-0.99 + 1.31 * gfc_stars[[1]][ind_vals[x],]) * t(p_ldd_mamm * mamm_range_mat)))
  cell_prod0_mamm <- colsums(preds_df_mamm$p_ldd0_mamm * mamm_range_mat)



  bird_range_mat <- st_intersects(bird_ranges_sf, rc_sf) %>% as.matrix()

  # get the land-use dependent probability of long distance dispersal at each cell
  p_ldd_bird <- 1-plnorm(q = ldd_dist_bird,
                         meanlog = log(preds_df_bird$bird_meanlog_inside %*%
                                         t(exp(-1.86 * fp_stars[[1]][ind_vals[x],]/50))),
                         sdlog = my_sdlog)

  # Get the product
  cell_prod_bird <- colsums(t(plogis(-0.99 + 1.31 * gfc_stars[[1]][ind_vals[x],]) * t(p_ldd_bird * bird_range_mat)))
  cell_prod0_bird <- colsums(preds_df_bird$p_ldd0_bird * bird_range_mat)

  (cell_prod_mamm + cell_prod_bird) / (cell_prod0_mamm + cell_prod0_bird)

})
time_3 <- Sys.time()
time_3
time_3 - time_2


# Save this file
setwd(output_dir)
saveRDS(out_list, file = output_filename)


"),fill=TRUE)
  sink()

  submit_script_name <- paste0("submit.script.2020.", i)
  sink(submit_script_name)
  cat(paste0("#!/bin/bash
#SBATCH -n 16 #request cpus
#SBATCH -N 1 #request nodes
#SBATCH -t 0-08:00 #time
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
