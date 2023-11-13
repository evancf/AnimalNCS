## AnimalNCS
### Assessing the role of animals in natural climate solutions

The goal of this analysis is to understand how disruption of animal-mediated seed dispersal influences the carbon accumulation potential of regenerating forests in the tropics. 


#### Explanation of scripts
The following describes the overall goals of each script in the analysis folder, which are intended to be run in their numeric order:

* 01_spatial_data_from_gee.R - This script downloads various spatial data to be used in the analyses. Before executing this script, users need to get set up to use Google Earth Engine in R via the 'rgee' package; see setup notes below. The script gives details on additional data downloads necessary here; these are also summarized in individual "Readme.txt" files within the sub-folders of the "./data/spatial data/" folder.
* 02_biodiv_data_prep.R - This script pulls together data from Cook-Patton et al. 2019 and Bukoski et al. 2022 on forest carbon accumulation and obtains values for covariates related to forest growth at these sites.
* 03_get_species_presence_by_sitecoords.R - This script determines which bird and mammal species are in range at each of the forest regeneration sites. Note that users must request range data from IUCN Red List / Bird Life International in order to reproduce this script.
* 04_land_use_dependent_dispersal_scripts - Note that this is a folder with an additional set of scripts contained within it. The goals of these scripts are to evaluate how human activities influence animal presence and movement and thus seed dispersal function.
    * 00_movebank_login.R - This script simply logs the user in to their Movebank account. Note that users must have an account with Movebank in order to execute this and the next scripts in this folder. Registration information here: https://www.movebank.org/cms/movebank-main.
    * 01_download_data.R - This script downloads available data on animal movement from Movebank.
    * 02_data_manipulation - This script subsets and reorganizes data to create a dataset of movement tracks over the course of a typical seed dispersal event.
    * 03_munge_add_covariates.R - This script adds covariates to the reorganized data and performs miscellaneous cleaning.
    * 03a_gee_covariates.R - This script adds covariates from the Earth Engine data catalog.
    * 03b_add_spatial_covariates.R - This script adds further spatial covariates.
    * 04_displacement_model.R - This script analyzes the data on animal movement prepared by the above scripts within this folder.
    * 20_fragmentation_model.R - This script analyzes the effect of fragmentation on species occurrence.
    * 30_plotting_model_results.R - This script produces a figure summarizing the above results of the displacement and fragmentation models.
* 05_get_dispersal_estimates.R - This script estimates seed dispersal function at each forest regeneration site (whose site data is assembled in 02_biodiv_data_prep.R) given local species composition, the seed dispersal function each species provides, and landscape factors that cause seed dispersal disruption. Note that this script also calls the script "disp_estimates_Fricke2022.R" to perform calculations based on the results of Fricke et al. 2022, "The effects of defaunation on plants' capacity to track climate change".
* 06_biodiv_analysis.R - This script presents the analysis of the effect of seed dispersal disruption on carbon accumulation. An equivalent model is fitted without the seed dispersal variable.
* 07_dispersal_mode_figure.R - This script develops a dataset on seed dispersal mode, estimates the prevalence of biotic dispersal across vegetation plots in a global database of vegetation plots, and then fits a random forest model to map the prevalence of biotic seed dispersal across forest and savanna biomes.
* 08_calculate_dispersal_integrity_on_cluster.R - This script outputs scripts that enable parallel computing of seed dispersal integrity (essentially the opposite of seed dispersal disruption) on a cluster. When executed on a cluster, these scripts produces tiles of seed dispersal integrity that together make a global map.
* 09_natural_regrowth_estimates.R - This script uses results of the fitted model from 06_biodiv_analysis.R to estimate natural regrowth rates under alternate scenarios (e.g., with or without seed dispersal disruption), in different areas (e.g., in restoration sites identified by Griscom et al. 2020 or across all tropical sites), and at two epochs (2000 and 2020).
* 10_natural_regrowth_estimates_uncertainty_on_cluster.R - This script produces scripts to be run in parallel on a HPC cluster. These sample from the posterior distribution of the fitted model (from 06_biodiv_analysis.R) and map natural regrowth 100 times, each time producing a map of regrowth estimates.
* 11_uncertainty_metrics_on_cluster.R - This script produces scripts to be run in parallel on a HPC cluster. This develops maps produced from 100 posterior samples of the fited model described above.
* 12_uncertainty_build_map.R - This script takes the maps produced by the script above and calculates credible intervals to show global maps of uncertainty.
* 13_natural_regrowth_metrics_plots.R - This script calculates summary metrics of carbon accumulation potential, their uncertainty, and produces the carbon accumulation maps that are presented in the study.




#### Additional software setup

##### Google Earth Engine.
Several scripts download data from the Earth Engine catalog using the 'rgee' package. Before being able to run these scripts, users will need to have signed up for a Google Earth Engine account and install rgee and its dependencies. Detailed information on getting Earth Engine access can be found here: https://developers.google.com/earth-engine/guides/access. Detailed information on getting set up to use Earth Engine via R using the rgee package is provided here: https://github.com/r-spatial/rgee. Note that this resource outlines further dependencies that need to be installed or set up in order to use the rgee package.

##### Spatial R packages
The scripts use R packages including 'sf' and 'terra' for spatial data. If a user has not used these packages before, follow the installation instructions here to install supporting software: https://r-spatial.github.io/sf/.

##### JAGS
For the Bayesian analysis, the scripts use JAGS (Just Another Gibbs Sampler), which needs to be installed before being able to use it via the 'rjags' package in R. Download information can be found here: https://mcmc-jags.sourceforge.io/.


#### Data notes

All datasets used in the analyses are publicly available, and the sources are described in the Data Availability section of the manuscript. Wherever possible, the scripts contain code to download the data directly from their sources. There are cases where this is not possible because the data can only be distributed to users after an individual submits a request to the data owner or because data are not shared in a way that enables downloading with code. In these cases, the data are either provided in this data and code package or comments within the script point users to the data source for individual download. The Readme.txt files within sub-folders of the data directory give more details regarding the locations where certain files should be placed.


#### Computing notes

Most of the analyses were executed on an Apple Macbook Pro M2 Max with 32 GB memory. Comments within the scripts sometimes offer estimates of the time required to execute portions of the scripts, but please note that this will vary depending on the user's computer. Data download speeds are also likely to depend on server load.

The most computationally intensive parts of the analysis were executed on a compute cluster, as noted within the scripts. The scripts with the suffix "on_cluster.R" themselves write out multiple R scripts and associated job scripts that can be used to execute pieces of the computation in parallel. Users can run these scripts on a typical university HPC cluster.
