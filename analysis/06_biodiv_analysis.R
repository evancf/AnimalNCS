source(list.files("./R", full.names = T))

ipak(c("tidyverse"))

# Joining datasets -------------------------------------------------------------

all_dat <- read.csv("./data/all_dat.csv", row.names = 1) %>% tibble()
coord_covars <- read.csv("./data/coord_covars.csv", row.names = 1) %>% tibble()
disp_covars <- read.csv("./data/disp_covars.csv", row.names = 1) %>% tibble()
ecoreg_covars <- read.csv("./data/ecoreg_dat.csv", row.names = 1) %>% tibble()

c_dat <- all_dat %>%
  filter(monoculture == "yes" | variables_name %in% c("aboveground_carbon",
                                                      "aboveground_biomass")) %>%
  mutate(variables_name = ifelse(monoculture == "yes", "aboveground_carbon", variables_name)) %>%
  mutate(mean_ha = ifelse(grepl("biomass", variables_name), mean_ha * 0.47, mean_ha)) %>%
  left_join(coord_covars) %>%
  left_join(disp_covars) %>%
  left_join(ecoreg_covars %>%
              dplyr::select(coord_id, BIOME_NAME)) %>%
  mutate(ratio = disp_both_eff_le_sum /disp_both_eff_no_le_sum) %>%
  filter(stand_age > 0, # Because we're log-transforming, these time = 0 will be removed
        mean_ha > 0)

# Ensuring the biome is the same across the Bukoski and our own attribution
biome_test <- c_dat %>% filter(!is.na(biome)) %>% dplyr::select(biome, BIOME_NAME)
all(biome_test[,1] == biome_test[,2])

# Remove the duplicate column
c_dat <- c_dat %>%
  dplyr::select(-biome) %>%
  mutate(biome = BIOME_NAME) %>%
  dplyr::select(-BIOME_NAME)
table(is.na(c_dat$biome))

# Filter to tropical biomes
c_dat <- c_dat %>%
  filter(grepl("ropical", biome))


# Get data object together for jags --------------------------------------------

j_dat <- c_dat %>%
  dplyr::select(site_id,
                plot_id,
                mean_ha,
                stand_age,
                ratio,
                monoculture,
                tap_mm,
                mat_C,
                npp_10000_max_average,
                lat_dec) %>%
  filter(complete.cases(.)) %>%
  mutate(DISP = ratio, # This is the dispersal integrity value
         MONO = ifelse(monoculture == "yes", 1, 0),
         NPP = npp_10000_max_average,
         LAT = abs(lat_dec)/10,
         AMP = tap_mm,
         AMT = mat_C,
         SITE = as.numeric(factor(site_id)))

n_site <- length(unique(j_dat$SITE))

# It's critical that dataframe is ordered by site
j_dat <- j_dat[order(j_dat$SITE),]

# Here is the data object given to jags
dataList <- list(
  BIOMASS = j_dat$mean_ha,
  AGE = j_dat$stand_age,

  n = length(j_dat$mean_ha),
  n_site = n_site,

  SITE = j_dat$SITE,

  MONO_for_SD = (j_dat$monoculture %>% as.factor() %>% as.numeric()),

  # DISP
  c1 = j_dat[!duplicated(j_dat$SITE),]$DISP,
  # AMP
  c2 = j_dat[!duplicated(j_dat$SITE),]$AMP %>% scale() %>% as.vector(),
  # AMT
  c3 = j_dat[!duplicated(j_dat$SITE),]$AMT %>% scale() %>% as.vector(),
  # NPP
  c4 = j_dat[!duplicated(j_dat$SITE),]$NPP %>% scale() %>% as.vector(),

  MONO = j_dat[!duplicated(j_dat$SITE),]$MONO %>% as.vector(),
  MONO_factor = j_dat[!duplicated(j_dat$SITE),]$monoculture %>% as.factor() %>% as.numeric(),

  pred_ages = 1:40,
  n_pred_ages = length(1:40),
  pred_mono_levels = c(0, 1),
  pred_disp_levels = c(0, 0.58),

  pred_disp_values = seq(0, 0.58, length.out = 30),
  n_pred_disp_values = 30
)


cov_to_var <- c("DISP", "AMP", "AMT", "NPP")
cov_to_name <- c("DISP", "AMP", "AMT", "NPP")
mono_intx <- c(1) # Which of the variables should have interactions with monoculture?

cov_to_sample <- c("c1", "c2", "c3", "c4")

n_cov <- length(cov_to_sample)

a_vars_to_sample <- paste("beta_a_",  cov_to_sample, sep = "")
b_vars_to_sample <- paste("beta_b_",  cov_to_sample, sep = "")

if(length(mono_intx) > 0){
  a_vars_to_sample <- c(a_vars_to_sample, "beta_a_mono", paste0("beta_a_mono_intx_c", mono_intx))
  b_vars_to_sample <- c(b_vars_to_sample, "beta_b_mono", paste0("beta_b_mono_intx_c", mono_intx))
}


# # The chunk below takes hours to run, so will just read in
# # a saved version of model outputs below.
# # Define model
#
# sink("./analysis/biomass_jags.txt")
# cat(paste0("
#
# model {
#
#   # Likelihood
#
#   for(i in 1:n){
#
#     # Biomass values lognormally distrubuted around a Monod function
#     BIOMASS[i] ~ dlnorm(log(mean[i]), tau_log_bio[MONO_for_SD[i]])
#     mean[i] <- a[SITE[i]] * AGE[i] / (b[SITE[i]] + AGE[i])
#
#   }
#
#   for(i in 1:n_site){
#     a[i] <- exp(beta_a_0 + beta_a_c1 * c1[i] +  beta_a_c2 * c2[i] +  beta_a_c3 * c3[i] +  beta_a_c4 * c4[i] + beta_a_mono * MONO[i] + beta_a_mono_intx_c1 * MONO[i] * c1[i] +
#
#       rand_a_s[SITE[i]])
#
#     b[i] <- exp(beta_b_0 + beta_b_c1 * c1[i] +  beta_b_c2 * c2[i] +  beta_b_c3 * c3[i] +  beta_b_c4 * c4[i] + beta_b_mono * MONO[i] + beta_b_mono_intx_c1 * MONO[i] * c1[i] +
#
#       rand_b_s[SITE[i]])
#
#   }
#
#
#   # Priors
#
#   beta_a_0 ~ dnorm(5, 0.1) # Roughly reasonable max carbon estimates
#
#   beta_a_c1 ~ dnorm(0, 0.1)
#   beta_a_c2 ~ dnorm(0, 0.1)
#   beta_a_c3 ~ dnorm(0, 0.1)
#   beta_a_c4 ~ dnorm(0, 0.1)
#
#   beta_a_mono ~ dnorm(0, 0.1)
#   beta_a_mono_intx_c1 ~ dnorm(0, 0.1)
#
#   beta_b_0 ~ dnorm(3, 0.1) # Using reasonable times to half saturation
#
#   beta_b_c1 ~ dnorm(0, 0.1)
#   beta_b_c2 ~ dnorm(0, 0.1)
#   beta_b_c3 ~ dnorm(0, 0.1)
#   beta_b_c4 ~ dnorm(0, 0.1)
#
#   beta_b_mono ~ dnorm(0, 0.1)
#   beta_b_mono_intx_c1 ~ dnorm(0, 0.1)
#
#
#
#   # Random effect for prior use, with sum to zero constraints
#   for(i in 1:(n_site-1)){
#     rand_a_s[i] ~ dnorm(0, tau_a_s)
#     rand_b_s[i] ~ dnorm(0, tau_b_s)
#   }
#   rand_a_s[n_site] <- -sum(rand_a_s[1:(n_site-1)])
#   rand_b_s[n_site] <- -sum(rand_b_s[1:(n_site-1)])
#
#   tau_a_s <- pow(sig_a_s, -2)
#   sig_a_s ~ dunif(0,10)
#   tau_b_s <- pow(sig_b_s, -2)
#   sig_b_s ~ dunif(0,10)
#
#   for(i in 1:2){
#     tau_log_bio[i] <- pow(sig_log_bio[i], -2)
#     sig_log_bio[i] ~ dunif(0, 10)
#   }
#
#
#
#   # Derived quantities
#
#   for(j in 1:2){ # Levels of dispersal to test
#     for(k in 1:2){ # Levels of monculture/non-monoculture
#
#       # Note that we are setting all values for other predictor variables to the
#       # mean value (because centered and scaled, this is zero)
#
#       aa[j,k] <- exp(beta_a_0 +  beta_a_c1 * pred_disp_levels[j] + beta_a_mono * pred_mono_levels[k] + beta_a_mono_intx_c1 * pred_mono_levels[k] * pred_disp_levels[j])
#
#       bb[j,k] <- exp(beta_b_0 +  beta_b_c1 * pred_disp_levels[j] + beta_b_mono * pred_mono_levels[k] + beta_b_mono_intx_c1 * pred_mono_levels[k] * pred_disp_levels[j])
#
#       for(t in 1:n_pred_ages){
#         pred_mean[j,k,t] <- aa[j,k] * pred_ages[t] / (bb[j,k] + pred_ages[t])
#
#         pred_biomass[j,k,t] ~ dlnorm(log(pred_mean[j,k,t]), tau_log_bio[k])
#
#       }
#     }
#   }
#
#   # Estimates at 30 years of growth
#   for(j in 1:n_pred_disp_values){ # Values of dispersal to test
#     for(k in 1:2){ # Levels of monculture/non-monoculture
#
#       # Note that we are setting all values for other predictor variables to the
#       # mean value (because centered and scaled, this is zero)
#
#       aa30[j,k] <- exp(beta_a_0 +  beta_a_c1 * pred_disp_values[j] + beta_a_mono * pred_mono_levels[k] + beta_a_mono_intx_c1 * pred_mono_levels[k] * pred_disp_values[j])
#
#       bb30[j,k] <- exp(beta_b_0 +  beta_b_c1 * pred_disp_values[j] + beta_b_mono * pred_mono_levels[k] + beta_b_mono_intx_c1 * pred_mono_levels[k] * pred_disp_values[j])
#
#       pred_mean30[j,k] <- aa30[j,k] * 30 / (bb30[j,k] + 30) # Monod
#
#     }
#   }
#
#
# } # End of model
#
# "),fill=TRUE)
# sink()
#
#
# # Will comment this out and pull in the saved output below because it
# # takes a while to run and requires jags installation
# ipak("rjags")
#
# set.seed(4)
# biomass_jags <- jags.model("analysis/biomass_jags.txt", data = dataList,
#                            n.adapt = 10000,
#                            n.chains = 3)
#
# update(biomass_jags, 10000)
#
# variables_to_sample <- c("beta_a_0",
#                          a_vars_to_sample,
#
#                          "beta_b_0",
#                          b_vars_to_sample,
#
#                          "rand_a_s",
#                          "rand_b_s",
#                          "pred_mean",
#                          "pred_biomass",
#
#                          "pred_mean30",
#                          "sig_log_bio"
# )
#
# biomass_samples <- coda.samples(biomass_jags,
#                                 variable.names = variables_to_sample,
#                                 thin = 500,
#                                 n.iter = 500000)
#
# # Save biomass samples and a few other helper vectors
# # Will include attributes for centering/scaling
# scale_df <- tibble(coef = "amt",
#                    center = attributes(scale(j_dat[!duplicated(j_dat$SITE),]$AMT))[[2]],
#                    scale = attributes(scale(j_dat[!duplicated(j_dat$SITE),]$AMT))[[3]]) %>%
#   bind_rows(tibble(coef = "tap",
#                    center = attributes(scale(j_dat[!duplicated(j_dat$SITE),]$AMP))[[2]],
#                    scale = attributes(scale(j_dat[!duplicated(j_dat$SITE),]$AMP))[[3]])) %>%
#   bind_rows(tibble(coef = "npp",
#                    center = attributes(scale(j_dat[!duplicated(j_dat$SITE),]$NPP))[[2]],
#                    scale = attributes(scale(j_dat[!duplicated(j_dat$SITE),]$NPP))[[3]]))
#
# # Actually will be easier if the rownames are the var_names
# scale_df <- as.data.frame(scale_df)
# rownames(scale_df) <- scale_df$coef
# scale_df <- scale_df[,-1]
#
# save(file = "./outputs/biomass_samples.RData",
#      biomass_samples,
#      variables_to_sample,
#      cov_to_var,
#      cov_to_name,
#      cov_to_sample,
#      n_cov,
#      a_vars_to_sample,
#      b_vars_to_sample,
#      scale_df)
#
#
# # Examine chains
# MCMCvis::MCMCtrace(biomass_samples,
#                    params = c(
#                      "sig_log_bio",
#
#                      "beta_a_0",
#                      a_vars_to_sample,
#
#                      "beta_b_0",
#                      b_vars_to_sample#,
#                      #"rand_a_s",
#                      #"rand_b_s"
#                    ),
#                    ISB = T,
#                    Rhat = T,
#                    n.eff = T,
#                    ind = T,
#                    pdf = T)
#
# MCMCvis::MCMCplot(biomass_samples,
#                   params = c("beta_a_0",
#                              a_vars_to_sample,
#                              "beta_b_0",
#                              b_vars_to_sample),
#         ci = c(50, 95))



# # Conduct a similar analysis but without the dispersal variable in the model.
# # This serves to fit a model of forest growth that uses typical environmental
# # variables.
#
# # Again, will save the model output and comment this out
#
# nd_cov_to_var <- c("AMP", "AMT", "NPP")
# nd_cov_to_name <- c("AMP", "AMT", "NPP")
#
# nd_cov_to_sample <- c("c2", "c3", "c4")
#
# n_nd_cov <- length(nd_cov_to_sample)
#
# a_vars_to_sample_nd <- paste("beta_a_",  nd_cov_to_sample, sep = "")
# b_vars_to_sample_nd <- paste("beta_b_",  nd_cov_to_sample, sep = "")
#
#
# # Define model
#
# sink("./analysis/biomass_nd_jags.txt")
# cat(paste0("
#
# model {
#
#   # Likelihood
#
#   for(i in 1:n){
#
#     # Biomass values lognormally distrubuted around a Monod function
#     BIOMASS[i] ~ dlnorm(log(mean[i]), tau_log_bio[MONO_for_SD[i]])
#     mean[i] <- a[SITE[i]] * AGE[i] / (b[SITE[i]] + AGE[i])
#
#   }
#
#   for(i in 1:n_site){
#     a[i] <- exp(beta_a_0 + beta_a_c2 * c2[i] +  beta_a_c3 * c3[i] +  beta_a_c4 * c4[i] + beta_a_mono * MONO[i] +
#
#       rand_a_s[SITE[i]])
#
#     b[i] <- exp(beta_b_0 + beta_b_c2 * c2[i] +  beta_b_c3 * c3[i] +  beta_b_c4 * c4[i] + beta_b_mono * MONO[i] +
#
#       rand_b_s[SITE[i]])
#
#   }
#
#
#   # Priors
#
#   beta_a_0 ~ dnorm(5, 0.1) # Roughly reasonable max carbon estimates
#
#   #beta_a_c1 ~ dnorm(0, 0.1)
#   beta_a_c2 ~ dnorm(0, 0.1)
#   beta_a_c3 ~ dnorm(0, 0.1)
#   beta_a_c4 ~ dnorm(0, 0.1)
#
#   beta_a_mono ~ dnorm(0, 0.1)
#   #beta_a_mono_intx_c1 ~ dnorm(0, 0.1)
#
#   beta_b_0 ~ dnorm(3, 0.1) # Using reasonable times to half saturation
#
#   #beta_b_c1 ~ dnorm(0, 0.1)
#   beta_b_c2 ~ dnorm(0, 0.1)
#   beta_b_c3 ~ dnorm(0, 0.1)
#   beta_b_c4 ~ dnorm(0, 0.1)
#
#   beta_b_mono ~ dnorm(0, 0.1)
#   #beta_b_mono_intx_c1 ~ dnorm(0, 0.1)
#
#
#
#   # Random effect for prior use, with sum to zero constraints
#   for(i in 1:(n_site-1)){
#     rand_a_s[i] ~ dnorm(0, tau_a_s)
#     rand_b_s[i] ~ dnorm(0, tau_b_s)
#   }
#   rand_a_s[n_site] <- -sum(rand_a_s[1:(n_site-1)])
#   rand_b_s[n_site] <- -sum(rand_b_s[1:(n_site-1)])
#
#   tau_a_s <- pow(sig_a_s, -2)
#   sig_a_s ~ dunif(0,10)
#   tau_b_s <- pow(sig_b_s, -2)
#   sig_b_s ~ dunif(0,10)
#
#   for(i in 1:2){
#     tau_log_bio[i] <- pow(sig_log_bio[i], -2)
#     sig_log_bio[i] ~ dunif(0, 10)
#   }
#
# } # End of model
#
# "),fill=TRUE)
# sink()
#
#
# # Will comment this out and pull in the saved output below because it
# # takes a while to run and requires jags installation
# ipak("rjags")
#
# set.seed(4)
# biomass_nd_jags <- jags.model("analysis/biomass_nd_jags.txt", data = dataList,
#                               n.adapt = 10000,
#                               n.chains = 3)
#
# update(biomass_nd_jags, 10000)
#
# variables_to_sample_nd <- c("beta_a_0",
#                             a_vars_to_sample_nd,
#
#                             "beta_b_0",
#                             b_vars_to_sample_nd,
#
#                             "rand_a_s",
#                             "rand_b_s",
#
#                             "sig_log_bio"
# )
#
# biomass_nd_samples <- coda.samples(biomass_nd_jags,
#                                    variable.names = variables_to_sample_nd,
#                                    thin = 500,
#                                    n.iter = 500000)
#
# # Save biomass samples and a few other helper vectors
# save(file = "./outputs/biomass_nd_samples.RData",
#      biomass_nd_samples,
#      variables_to_sample_nd,
#      nd_cov_to_var,
#      nd_cov_to_name,
#      nd_cov_to_sample,
#      n_nd_cov,
#      a_vars_to_sample_nd,
#      b_vars_to_sample_nd)
#
#
# # Examine chains
# MCMCvis::MCMCtrace(biomass_nd_samples,
#                    params = c(
#                      "sig_log_bio",
#
#                      "beta_a_0",
#                      a_vars_to_sample_nd,
#
#                      "beta_b_0",
#                      b_vars_to_sample_nd#,
#                      #"rand_a_s",
#                      #"rand_b_s"
#                    ),
#                    ISB = T,
#                    Rhat = T,
#                    n.eff = T,
#                    ind = T,
#                    pdf = T)






# Figures ----------------------------------------------------------------------

# Load posterior samples if not run above
load(file = "./outputs/biomass_samples.RData")

# Get posterior means and credible interval vals
biomass_samples_df <- do.call(rbind, biomass_samples)
bsm <- colMeans(biomass_samples_df)
bs975 <-apply(biomass_samples_df, 2, function(x) quantile(x, 0.975))
bs025 <-apply(biomass_samples_df, 2, function(x) quantile(x, 0.025))

# Get the names of these coefs
beta_vec <- names(bsm)[grepl("beta", names(bsm))]
pred_vec <- names(bsm)[grepl("pred_biomass", names(bsm))]
mean_vec <- names(bsm)[grepl("pred_mean", names(bsm))]




# A couple functions to help plotting the results
sample_posterior <- function(chains){
  cc <- sample(1:length(chains), 1)
  ii <- sample(1:length(chains[[1]][,1]), 1)

  chains[[cc]][ii, beta_vec]
}

curve_intx <- function(mean_or_sample = "mean",
                       x,
                       c1_val = 0,
                       c2_val = 0,
                       c3_val = 0,
                       c4_val = 0,
                       MONO_val = 0,
                       ...){

  if(mean_or_sample == "sample"){
    x <- sample_posterior(biomass_samples) %>% t() %>% as.data.frame()
  }

  a_val <- exp(x$beta_a_0 + x$beta_a_c1 * c1_val +  x$beta_a_c2 * c2_val +
                 x$beta_a_c3 * c3_val + x$beta_a_c4 * c4_val +
                 x$beta_a_mono * MONO_val +
                 x$beta_a_mono_intx_c1 * MONO_val * c1_val)

  b_val <- exp(x$beta_b_0 + x$beta_b_c1 * c1_val +  x$beta_b_c2 * c2_val +
                 x$beta_b_c3 * c3_val + x$beta_b_c4 * c4_val +
                 x$beta_b_mono * MONO_val +
                 x$beta_b_mono_intx_c1 * MONO_val * c1_val)

  curve(
    a_val * x / (b_val + x),
    ylim = c(0, 200),
    xlim = c(0, 40),
    ...)

}



png(file = "./outputs/figures/biomass model results.png",
    units = "in",
    width = 4.76,#7.24,
    height = 2.6,#6,
    res = 440)

par(mfrow=c(1,2))
par(mar = c(3.5, 3.9, 0.1, 0.1))
set.seed(111)

# Plot aboveground biomass over stand age based on model predictions that
# account for site-level non-independence and error

j_dat$pred_a <- exp(bsm["beta_a_0"] +
                      bsm["beta_a_c1"] * j_dat$DISP +
                      bsm["beta_a_c2"] * ((j_dat$AMP - scale_df["tap", "center"]) / scale_df["tap", "scale"]) +
                      bsm["beta_a_c3"] * ((j_dat$AMT - scale_df["amt", "center"]) / scale_df["amt", "scale"]) +
                      bsm["beta_a_c4"] * ((j_dat$NPP - scale_df["npp", "center"]) / scale_df["npp", "scale"]))
j_dat$pred_b <- exp(bsm["beta_b_0"] +
                      bsm["beta_b_c1"] * j_dat$DISP +
                      bsm["beta_b_c2"] * ((j_dat$AMP - scale_df["tap", "center"]) / scale_df["tap", "scale"]) +
                      bsm["beta_b_c3"] * ((j_dat$AMT - scale_df["amt", "center"]) / scale_df["amt", "scale"]) +
                      bsm["beta_b_c4"] * ((j_dat$NPP - scale_df["npp", "center"]) / scale_df["npp", "scale"]))
j_dat$pred_est <- j_dat$pred_a * j_dat$stand_age / (j_dat$pred_b + j_dat$stand_age)


legend_cex <- 0.8

polygon_alpha <- 0.8

mgp_vals_y <- c(3, 0.5, 0)
mgp_vals_x <- c(3, 0.4, 0)
tck_val <- -0.04

axis_cex <- 0.9

line_x <- 1.5
line_y <- 2

plot(NA,
     xlim = c(0, 50),
     ylim = c(0, 120),
     las = 1,
     xlab = "",#"Stand age (years)",
     ylab = "",#"Aboveground carbon (MgC/ha)",
     xaxt = "n",
     yaxt = "n",
     frame = F,
     mgp = mgp_vals_y,
     tck = tck_val)
mtext("Stand age", side = 1, line = line_x)
mtext("Aboveground\ncarbon (MgC/ha)", side = 2, line = line_y)

text(-25.5, 120, "a", xpd = T,
     font = 2)

axis(1, at = c(0, 25, 50),
     labels = c("0", "25", ">50"),
     mgp = mgp_vals_x,
     tck = tck_val,
     cex.axis = axis_cex)
axis(2, at = c(0, 40, 80, 120),
     labels = c("0", "40", "80", ">120"),
     mgp = mgp_vals_y,
     tck = tck_val,
     las = 1,
     cex.axis = axis_cex)

color_range <- colorRampPalette(c("red", "purple", "blue"))(100)
j_plot <- j_dat %>% #j_dat[sample(1:nrow(j_dat)),] %>%
  arrange(ratio) %>%
  filter(monoculture == "no") %>%
  filter(!duplicated(paste(site_id, stand_age)))
points(x = ifelse(j_plot$stand_age > 50, 50, j_plot$stand_age),
       y = ifelse(j_plot$pred_est > 120, 120, j_plot$pred_est),
  col = alpha(color_range[findInterval(j_plot$ratio, seq(0, 0.58, length.out = 100))], 0.5),
  pch = 16,
  cex = 0.5)

# A super clunky way to show a gradient legend
library("fields")
par(new = TRUE, pty = "m", pin = c(0.2, 1))
image.plot(zlim = c(-1, 1),
           col = alpha(rev(color_range), 0.7),
           horizontal = T,
           legend.only = TRUE,
           legend.width = 0.3,
           legend.shrink = 3,
           axis.args = list(xaxt = "n"),
           legend.mar = 11)
text(25, 231,
     pos = 1,
     "Dispersal\ndisruption",
     col = "grey30",
     cex = legend_cex,
     font = 3)
text(x = c(-55, 25, 105),
     y = 166,
     labels = c("0", "0.5", "1"),
     cex = legend_cex,
     col = "grey30")




# # Here's some other code that just looks at the ends of the spectrum.
# # Decided this wasn't as effective a representation of the model, so used
# # the above code instead.
# lo_disp_col <- rgb(239,138,98, maxColorValue = 255)
# hi_disp_col <- rgb(103,169,207, maxColorValue = 255)
# legend_cex <- 0.8
#
# polygon_alpha <- 0.8
#
# mgp_vals_y <- c(3, 0.5, 0)
# mgp_vals_x <- c(3, 0.4, 0)
# tck_val <- -0.04
#
# line_x <- 1.5
# line_y <- 2
#
# plot(NA,
#      xlim = c(0, 40),
#      ylim = c(0, 180),
#      las = 1,
#      xlab = "",#"Stand age (years)",
#      ylab = "",#"Aboveground carbon (MgC/ha)",
#      xaxt = "n",
#      yaxt = "n",
#      frame = F,
#      mgp = mgp_vals_y,
#      tck = tck_val)
# mtext("Stand age", side = 1, line = line_x)
# mtext("Aboveground carbon\n(MgC/ha)", side = 2, line = line_y)
#
# axis(1, at = c(0, 20, 40),
#      mgp = mgp_vals_x,
#      tck = tck_val)
# axis(2, at = c(0, 60, 120, 180),
#      mgp = mgp_vals_y,
#      tck = tck_val,
#      las = 1)
#
# polygon(y = c(0, bs975[grep("pred_mean[2,1,",
#                             names(bs975), fixed = T)]) %>%
#           c(rev(bs025[grep("pred_mean[2,1,",
#                            names(bs025), fixed = T)]), 0),
#         x = c(0, dataList$pred_ages) %>%
#           c(rev(c(0, dataList$pred_ages))),
#         col = alpha(hi_disp_col, polygon_alpha),
#         border = NA)
#
# lines(y = c(0, bsm[grep("pred_mean[2,1,",
#                         names(bsm), fixed = T)]),
#       x = c(0, dataList$pred_ages),
#       col = hi_disp_col,
#       lwd = 3,
#       lend = "butt")
#
# polygon(y = c(0, bs975[grep("pred_mean[1,1,",
#                             names(bs975), fixed = T)]) %>%
#           c(rev(bs025[grep("pred_mean[1,1,",
#                            names(bs025), fixed = T)]), 0),
#         x = c(0, dataList$pred_ages) %>%
#           c(rev(c(0, dataList$pred_ages))),
#         col = alpha(lo_disp_col, polygon_alpha),
#         border = NA)
#
# lines(y = c(0, bsm[grep("pred_mean[1,1,",
#                         names(bsm), fixed = T)]),
#       x = c(0, dataList$pred_ages),
#       col = lo_disp_col,
#       lwd = 3,
#       lend = "butt")
#
#
#
# # Legend
#
# text("Low\ndisruption\nscenario",
#      x = 31,
#      y = 155,
#      col = hi_disp_col,
#      font = 3,
#      cex = legend_cex,
#      pos = 2)
#
# text("High\ndisruption\nscenario",
#      x = 41,
#      y = 52,
#      col = lo_disp_col,
#      font = 3,
#      cex = legend_cex,
#      pos = 2)


#par(mar = c(5.1, 3.1, 4.1, 2.1))

# We want to put these in terms of dispersal disruption on a relative
# scale rather than dispersal integrity
dispersal_disruption_values <- (0.58 - dataList$pred_disp_values)/0.58

plot(NA,
     xlim = c(0, 1),
     ylim = c(0, 5),
     las = 1,
     xlab = "",#"Dispersal disruption",
     ylab = "",#"Accumulation rate (MgC/ha/yr)",
     frame = F,
     xaxt = "n",
     yaxt = "n",
     mgp = mgp_vals_y,
     tck = tck_val)
axis(1, at = c(0, 0.5, 1),
     labels = c(0, 0.5, 1),
     mgp = mgp_vals_x,
     tck = tck_val,
     cex.axis = axis_cex)
axis(2, at = c(0:5),
     las = 2,
     mgp = mgp_vals_y,
     tck = tck_val,
     cex.axis = axis_cex)

mtext("Seed dispersal\ndisruption", side = 1, line = line_x*1.67)

mtext("Accumulation\nrate (MgC/ha/yr)", side = 2, line = line_y/1.57)

mono_col <- rgb(216,179,101, maxColorValue = 255)
nat_col <- rgb(90,180,172, maxColorValue = 255)


polygon(y = (c(0, bs975[grepl("pred_mean30",
                              names(bs975), fixed = T) &
                          grepl(",2]",
                                names(bs975), fixed = T)]) %>%
               c(rev(bs025[grepl("pred_mean30",
                                 names(bs025), fixed = T) &
                             grepl(",2]",
                                   names(bs025), fixed = T)]), 0))/30,
        x = c(1, dispersal_disruption_values) %>%
          c(rev(c(1, dispersal_disruption_values))),
        col = alpha(mono_col, polygon_alpha * 0.6),
        border = NA)


lines(y = c(bsm[grepl("pred_mean30",
                      names(bsm),
                      fixed = T) &
                  grepl(",2]",
                        names(bsm),
                        fixed = T)]) / 30,
      x = c(dispersal_disruption_values),
      col = alpha(mono_col, polygon_alpha * 0.6),
      lwd = 3,
      lty = 2,
      lend = "butt")

polygon(y = (c(0, bs975[grepl("pred_mean30",
                              names(bs975), fixed = T) &
                          grepl(",1]",
                                names(bs975), fixed = T)]) %>%
               c(rev(bs025[grepl("pred_mean30",
                                 names(bs025), fixed = T) &
                             grepl(",1]",
                                   names(bs025), fixed = T)]), 0))/30,
        x = c(1, dispersal_disruption_values) %>%
          c(rev(c(1, dispersal_disruption_values))),
        col = alpha(nat_col, polygon_alpha),
        border = NA)


lines(y = c(bsm[grepl("pred_mean30",
                      names(bsm),
                      fixed = T) &
                  grepl(",1]",
                        names(bsm),
                        fixed = T)]) / 30,
      x = c(dispersal_disruption_values),
      col = nat_col,
      lwd = 3,
      lend = "butt")



text("Monoculture\nplantation",
     x = 0,
     y = 0.5,
     col = mono_col,
     font = 3,
     cex = legend_cex,
     pos = 4)

text("Natural\nregrowth",
     x = 0.1,
     y = 4.25,
     col = nat_col,
     font = 3,
     cex = legend_cex,
     pos = 4)

text(-0.41, 5, "b", xpd = T,
     font = 2)


dev.off()




# Effect sizes for npp, temp, precip, dispersal integrity (natural),
# dispersal integrity (monoculture) in terms of their effect on annualized
# carbon accumulation

# Note that temp precip npp were scaled prior to analysis as the values,
# allowing the values to be centered around zero. The dispersal integrity
# values weren't, so to put all these effect sizes into comparable units,
# will use a 1 sd unit. Also make this negative to make visualizations
# in terms of seed dispersal disruption.

# Below is the (kinda clunky) way we'll calculate the effect size on
# annualized growth rate.

aa1 <- list(c(),c(),c(),c(),c())
names(aa1) <- c("pNPP", "AMT", "TAP", "DD (natural)", "DD (monoculture)")
bb1 <- aa1
aa2 <- aa1
bb2 <- aa1
effect_size <- aa1

mean_c1 <- j_dat[!duplicated(j_dat$SITE),]$DISP %>% mean()
val_c1 <- -1 * j_dat[!duplicated(j_dat$SITE),]$DISP %>% sd()
val_c2 <- 1
val_c3 <- 1
val_c4 <- 1
val_mono <- 1

biomass_samples_df <- biomass_samples_df %>% as.data.frame()

# pNPP
aa1[["pNPP"]] <- exp(biomass_samples_df$beta_a_0 +
                       biomass_samples_df$beta_a_c1 * mean_c1 +
                       biomass_samples_df$beta_a_c2 * 0 +
                       biomass_samples_df$beta_a_c3 * 0 +
                       biomass_samples_df$beta_a_c4 * 0 +
                       biomass_samples_df$beta_a_mono * 0 +
                       biomass_samples_df$beta_a_mono_intx_c1 * 0 * mean_c1)

bb1[["pNPP"]] <- exp(biomass_samples_df$beta_b_0 +
                       biomass_samples_df$beta_b_c1 * mean_c1 +
                       biomass_samples_df$beta_b_c2 * 0 +
                       biomass_samples_df$beta_b_c3 * 0 +
                       biomass_samples_df$beta_b_c4 * 0 +
                       biomass_samples_df$beta_b_mono * 0 +
                       biomass_samples_df$beta_b_mono_intx_c1 * 0 * mean_c1)

aa2[["pNPP"]] <- exp(biomass_samples_df$beta_a_0 +
                       biomass_samples_df$beta_a_c1 * mean_c1 +
                       biomass_samples_df$beta_a_c2 * 0 +
                       biomass_samples_df$beta_a_c3 * 0 +
                       biomass_samples_df$beta_a_c4 * val_c4 +
                       biomass_samples_df$beta_a_mono * 0 +
                       biomass_samples_df$beta_a_mono_intx_c1 * 0 * mean_c1)

bb2[["pNPP"]] <- exp(biomass_samples_df$beta_b_0 +
                       biomass_samples_df$beta_b_c1 * mean_c1 +
                       biomass_samples_df$beta_b_c2 * 0 +
                       biomass_samples_df$beta_b_c3 * 0 +
                       biomass_samples_df$beta_b_c4 * val_c4 +
                       biomass_samples_df$beta_b_mono * 0 +
                       biomass_samples_df$beta_b_mono_intx_c1 * 0 * mean_c1)

effect_size[["pNPP"]] <- (aa2[["pNPP"]] * 30 / (bb2[["pNPP"]] + 30) / 30) -
  (aa1[["pNPP"]] * 30 / (bb1[["pNPP"]] + 30) / 30)
#hist(effect_size[["pNPP"]])




# AMT
aa1[["AMT"]] <- exp(biomass_samples_df$beta_a_0 +
                      biomass_samples_df$beta_a_c1 * mean_c1 +
                      biomass_samples_df$beta_a_c2 * 0 +
                      biomass_samples_df$beta_a_c3 * 0 +
                      biomass_samples_df$beta_a_c4 * 0 +
                      biomass_samples_df$beta_a_mono * 0 +
                      biomass_samples_df$beta_a_mono_intx_c1 * 0 * mean_c1)

bb1[["AMT"]] <- exp(biomass_samples_df$beta_b_0 +
                      biomass_samples_df$beta_b_c1 * mean_c1 +
                      biomass_samples_df$beta_b_c2 * 0 +
                      biomass_samples_df$beta_b_c3 * 0 +
                      biomass_samples_df$beta_b_c4 * 0 +
                      biomass_samples_df$beta_b_mono * 0 +
                      biomass_samples_df$beta_b_mono_intx_c1 * 0 * mean_c1)

aa2[["AMT"]] <- exp(biomass_samples_df$beta_a_0 +
                      biomass_samples_df$beta_a_c1 * mean_c1 +
                      biomass_samples_df$beta_a_c2 * 0 +
                      biomass_samples_df$beta_a_c3 * val_c3 +
                      biomass_samples_df$beta_a_c4 * 0 +
                      biomass_samples_df$beta_a_mono * 0 +
                      biomass_samples_df$beta_a_mono_intx_c1 * 0 * mean_c1)

bb2[["AMT"]] <- exp(biomass_samples_df$beta_b_0 +
                      biomass_samples_df$beta_b_c1 * mean_c1 +
                      biomass_samples_df$beta_b_c2 * 0 +
                      biomass_samples_df$beta_b_c3 * val_c3 +
                      biomass_samples_df$beta_b_c4 * 0 +
                      biomass_samples_df$beta_b_mono * 0 +
                      biomass_samples_df$beta_b_mono_intx_c1 * 0 * mean_c1)

effect_size[["AMT"]] <- (aa2[["AMT"]] * 30 / (bb2[["AMT"]] + 30) / 30) -
  (aa1[["AMT"]] * 30 / (bb1[["AMT"]] + 30) / 30)
#hist(effect_size[["AMT"]])




# TAP
aa1[["TAP"]] <- exp(biomass_samples_df$beta_a_0 +
                      biomass_samples_df$beta_a_c1 * mean_c1 +
                      biomass_samples_df$beta_a_c2 * 0 +
                      biomass_samples_df$beta_a_c3 * 0 +
                      biomass_samples_df$beta_a_c4 * 0 +
                      biomass_samples_df$beta_a_mono * 0 +
                      biomass_samples_df$beta_a_mono_intx_c1 * 0 * mean_c1)

bb1[["TAP"]] <- exp(biomass_samples_df$beta_b_0 +
                      biomass_samples_df$beta_b_c1 * mean_c1 +
                      biomass_samples_df$beta_b_c2 * 0 +
                      biomass_samples_df$beta_b_c3 * 0 +
                      biomass_samples_df$beta_b_c4 * 0 +
                      biomass_samples_df$beta_b_mono * 0 +
                      biomass_samples_df$beta_b_mono_intx_c1 * 0 * mean_c1)

aa2[["TAP"]] <- exp(biomass_samples_df$beta_a_0 +
                      biomass_samples_df$beta_a_c1 * mean_c1 +
                      biomass_samples_df$beta_a_c2 * val_c2 +
                      biomass_samples_df$beta_a_c3 * 0 +
                      biomass_samples_df$beta_a_c4 * 0 +
                      biomass_samples_df$beta_a_mono * 0 +
                      biomass_samples_df$beta_a_mono_intx_c1 * 0 * mean_c1)

bb2[["TAP"]] <- exp(biomass_samples_df$beta_b_0 +
                      biomass_samples_df$beta_b_c1 * mean_c1 +
                      biomass_samples_df$beta_b_c2 * val_c2 +
                      biomass_samples_df$beta_b_c3 * 0 +
                      biomass_samples_df$beta_b_c4 * 0 +
                      biomass_samples_df$beta_b_mono * 0 +
                      biomass_samples_df$beta_b_mono_intx_c1 * 0 * mean_c1)

effect_size[["TAP"]] <- (aa2[["TAP"]] * 30 / (bb2[["TAP"]] + 30) / 30) -
  (aa1[["TAP"]] * 30 / (bb1[["TAP"]] + 30) / 30)
#hist(effect_size[["TAP"]])




# DD (natural)
aa1[["DD (natural)"]] <- exp(biomass_samples_df$beta_a_0 +
                               biomass_samples_df$beta_a_c1 * 0 +
                               biomass_samples_df$beta_a_c2 * 0 +
                               biomass_samples_df$beta_a_c3 * 0 +
                               biomass_samples_df$beta_a_c4 * 0 +
                               biomass_samples_df$beta_a_mono * 0 +
                               biomass_samples_df$beta_a_mono_intx_c1 * 0 * 0)

bb1[["DD (natural)"]] <- exp(biomass_samples_df$beta_b_0 +
                               biomass_samples_df$beta_b_c1 * 0 +
                               biomass_samples_df$beta_b_c2 * 0 +
                               biomass_samples_df$beta_b_c3 * 0 +
                               biomass_samples_df$beta_b_c4 * 0 +
                               biomass_samples_df$beta_b_mono * 0 +
                               biomass_samples_df$beta_b_mono_intx_c1 * 0 * 0)

aa2[["DD (natural)"]] <- exp(biomass_samples_df$beta_a_0 +
                               biomass_samples_df$beta_a_c1 * val_c1 +
                               biomass_samples_df$beta_a_c2 * 0 +
                               biomass_samples_df$beta_a_c3 * 0 +
                               biomass_samples_df$beta_a_c4 * 0 +
                               biomass_samples_df$beta_a_mono * 0 +
                               biomass_samples_df$beta_a_mono_intx_c1 * 0 * val_c1)

bb2[["DD (natural)"]] <- exp(biomass_samples_df$beta_b_0 +
                               biomass_samples_df$beta_b_c1 * val_c1 +
                               biomass_samples_df$beta_b_c2 * 0 +
                               biomass_samples_df$beta_b_c3 * 0 +
                               biomass_samples_df$beta_b_c4 * 0 +
                               biomass_samples_df$beta_b_mono * 0 +
                               biomass_samples_df$beta_b_mono_intx_c1 * 0 * val_c1)

effect_size[["DD (natural)"]] <- (aa2[["DD (natural)"]] * 30 / (bb2[["DD (natural)"]] + 30) / 30) -
  (aa1[["DD (natural)"]] * 30 / (bb1[["DD (natural)"]] + 30) / 30)
#hist(effect_size[["DD (natural)"]])




# DD (monoculture)
aa1[["DD (monoculture)"]] <- exp(biomass_samples_df$beta_a_0 +
                                   biomass_samples_df$beta_a_c1 * 0 +
                                   biomass_samples_df$beta_a_c2 * 0 +
                                   biomass_samples_df$beta_a_c3 * 0 +
                                   biomass_samples_df$beta_a_c4 * 0 +
                                   biomass_samples_df$beta_a_mono * val_mono +
                                   biomass_samples_df$beta_a_mono_intx_c1 * val_mono * 0)

bb1[["DD (monoculture)"]] <- exp(biomass_samples_df$beta_b_0 +
                                   biomass_samples_df$beta_b_c1 * 0 +
                                   biomass_samples_df$beta_b_c2 * 0 +
                                   biomass_samples_df$beta_b_c3 * 0 +
                                   biomass_samples_df$beta_b_c4 * 0 +
                                   biomass_samples_df$beta_b_mono * val_mono +
                                   biomass_samples_df$beta_b_mono_intx_c1 * val_mono * 0)

aa2[["DD (monoculture)"]] <- exp(biomass_samples_df$beta_a_0 +
                                   biomass_samples_df$beta_a_c1 * val_c1 +
                                   biomass_samples_df$beta_a_c2 * 0 +
                                   biomass_samples_df$beta_a_c3 * 0 +
                                   biomass_samples_df$beta_a_c4 * 0 +
                                   biomass_samples_df$beta_a_mono * val_mono +
                                   biomass_samples_df$beta_a_mono_intx_c1 * val_mono * val_c1)

bb2[["DD (monoculture)"]] <- exp(biomass_samples_df$beta_b_0 +
                                   biomass_samples_df$beta_b_c1 * val_c1 +
                                   biomass_samples_df$beta_b_c2 * 0 +
                                   biomass_samples_df$beta_b_c3 * 0 +
                                   biomass_samples_df$beta_b_c4 * 0 +
                                   biomass_samples_df$beta_b_mono * val_mono +
                                   biomass_samples_df$beta_b_mono_intx_c1 * val_mono * val_c1)

effect_size[["DD (monoculture)"]] <- (aa2[["DD (monoculture)"]] * 30 / (bb2[["DD (monoculture)"]] + 30) / 30) -
  (aa1[["DD (monoculture)"]] * 30 / (bb1[["DD (monoculture)"]] + 30) / 30)
#hist(effect_size[["DD (monoculture)"]])


# Now produce a catepillar plot
png(file = "./outputs/figures/biomass model effect sizes.png",
    units = "in",
    width = 6.5,
    height = 4,#6,
    res = 440)

par(mar = c(5.1, 4.1, 1.1, 2.1))

x_vals <- 1:5
widths <- c(5,3,1)

plot(NA,
     xlim = c(0.25, 5.75),
     ylim = c(-0.5, 0.75),
     frame = F,
     xlab = "",
     ylab = "Accumulation rate ffect size",
     xaxt = "n",
     las = 2)
abline(h = 0,
       lty = 2)
text(x = x_vals + 0.2, y = -0.5,
     #labels = names(effect_size),
     labels = c("Potential NPP", "Temperature", "Precipitation",
                "Dispersal disruption\n(natural regrowth)",
                "Dispersal disruption\n(monoculture)"),
     cex = 0.8,
     pos = 2,
     srt = 45,
     xpd = T)

q50_lo <- lapply(effect_size, function(x) quantile(x, 0.25)) %>% unlist()
q50_hi <- lapply(effect_size, function(x) quantile(x, 0.75)) %>% unlist()
segments(x0 = x_vals,
         x1 = x_vals,
         y0 = q50_lo,
         y1 = q50_hi,
         lwd = widths[1],
         lend = "butt")

q95_lo <- lapply(effect_size, function(x) quantile(x, 0.025)) %>% unlist()
q95_hi <- lapply(effect_size, function(x) quantile(x, 0.975)) %>% unlist()
segments(x0 = x_vals,
         x1 = x_vals,
         y0 = q95_lo,
         y1 = q95_hi,
         lwd = widths[2],
         lend = "butt")


q99_lo <- lapply(effect_size, function(x) quantile(x, 0.005)) %>% unlist()
q99_hi <- lapply(effect_size, function(x) quantile(x, 0.995)) %>% unlist()
segments(x0 = x_vals,
         x1 = x_vals,
         y0 = q99_lo,
         y1 = q99_hi,
         lwd = widths[3],
         lend = "butt")


med_vals <- lapply(effect_size, median)
points(x_vals, med_vals,
       pch = c(16, 16, 16, 16, 21),
       cex = 1.2,
       bg = "white")

dev.off()










png(file = "./outputs/figures/figure s6.png",
    units = "in",
    width = 5,#7.24,
    height = 5,#6,
    res = 440)

site_dat <- j_dat %>% filter(!duplicated(j_dat$site_id)) %>%
  mutate(monoculture = fct_recode(monoculture,
                                  "Monoculture" = "yes",
                                  "Natural regrowth" = "no"))

p1 <- ggplot(site_dat, aes(x = monoculture, y = AMP)) +
  geom_violin(bw = diff(range(site_dat$AMP))/bw_val) +
  theme_classic() +
  labs(x = "", y = "Annual precipitation (mm)")

p2 <- ggplot(site_dat, aes(x = monoculture, y = AMT)) +
  geom_violin(bw = diff(range(site_dat$AMT))/bw_val) +
  theme_classic() +
  labs(x = "", y = "Annual mean temperature (°C)")

p3 <- ggplot(site_dat, aes(x = monoculture, y = NPP)) +
  geom_violin(bw = diff(range(site_dat$NPP))/bw_val) +
  theme_classic() +
  labs(x = "", y = "Potential NPP Index")

p4 <- ggplot(site_dat, aes(x = monoculture, y = lat_dec)) +
  geom_violin(bw = diff(range(site_dat$lat_dec))/bw_val) +
  theme_classic() +
  labs(x = "", y = "Latitude")

library("gridExtra")
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)

dev.off()
