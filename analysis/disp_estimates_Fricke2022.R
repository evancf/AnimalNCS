# This uses dispersal estimates from Fricke et al. 2022 Science, with the
# default argument values reflecting how these were presented in that
# original analysis. We include alternate ways to calculate dispersal
# estimates here, as well as include the ability to incorporate an effect
# of landscape modification on the dispersal estimates


# Going to make a function to calculate germination probabilities
# based on the log odds ratios predicted by the germ.mod.
# We will assume a 'baseline' germination of 40% based on the average
# germination percentage in the gut passage germination studies.

pgerm <- function(yi, prob = 0.4){
  (exp(yi) * prob) / (1 + exp(yi) * prob - prob)
}



dispersal_estimates <- function(dist_threshold_mamm = 1000,
                                dist_threshold_bird = 1000,
                                include_landscape_effect = F,
                                sum_or_mean = "sum",
                                type = "effective_dispersal"){

  # Assign a threshold for long-distance dispersal.
  # We will use 1000 meters
  ldd_dist_mamm <- dist_threshold_mamm
  ldd_dist_bird <- dist_threshold_bird
  my_sdlog <- 0.6 # Based on a meta analysis

  # We are going to make an object to describe the dispersal done by each bird
  # and mammal
  x_bird <- tibble(binomial = spp_names_bird)
  x_mamm <- tibble(binomial = spp_names_mamm)

  if(type == "effective_dispersal"){

    # Index of how many seeds dispersed. But critically, this also has to be
    # multiplied by fragmentation-dependent probability of occurrence. Empirically,
    # we found this relationship to be this: P(presence) = plogis(-0.99 + 1.31 * treecover)
    # when treecover is on a 0-1 scale
    if(include_landscape_effect){
      x_bird$mu <- particip_preds_bird *
        exp(quant_preds_bird) *
        int_preds_bird *
        pgerm(germ_preds_bird) %*% t(plogis(-0.99 + 1.31 * coord_covars$treecover2000/100))

      # P(presence) = plogis(-0.99 + 1.31 * treecover)

      # Portion of those seeds that fall greater than the long distance threshold.
      # But critically, this also has to be influenced by the reduction in
      # movement distances as a result of human footprint
      x_bird$p.ldd <- 1-plnorm(q = ldd_dist_bird,
                               meanlog = log(exp(gpt_preds_bird) *
                                               exp(move_preds_bird) *
                                               1000 / 60 *
                                               exp(corr_preds_bird) %*%
                                               t(exp(-1.86 * coord_covars$fp/50))),
                               sdlog = my_sdlog)
    } else{
      x_bird$mu <- particip_preds_bird *
        exp(quant_preds_bird) *
        int_preds_bird *
        pgerm(germ_preds_bird)
      x_bird$p.ldd <- 1-plnorm(q = ldd_dist_bird,
                               meanlog = log(exp(gpt_preds_bird) *
                                               exp(move_preds_bird) *
                                               1000 / 60 *
                                               exp(corr_preds_bird)),
                               sdlog = my_sdlog)
    }

    # The product of these two is the long-distance dispersal index
    x_bird$prod <- x_bird$p.ldd * x_bird$mu

    # The same for mammals...
    if(include_landscape_effect){

      x_mamm$mu <- particip_preds_mamm *
        exp(quant_preds_mamm) *
        int_preds_mamm *
        pgerm(germ_preds_mamm) %*% t(plogis(-0.99 + 1.31 * coord_covars$treecover2000/100))

      x_mamm$p.ldd <- 1-plnorm(q = ldd_dist_mamm,
                               meanlog = log(exp(gpt_preds_mamm) *
                                               exp(move_preds_mamm) *
                                               1000 / 60 *
                                               exp(corr_preds_mamm) %*%
                                               t(exp(-1.86 * coord_covars$fp/50))),
                               sdlog = my_sdlog)
    } else{
      x_mamm$mu <- particip_preds_mamm *
        exp(quant_preds_mamm) *
        int_preds_mamm *
        pgerm(germ_preds_mamm)

      x_mamm$p.ldd <- 1-plnorm(q = ldd_dist_mamm,
                               meanlog = log(exp(gpt_preds_mamm) *
                                               exp(move_preds_mamm) *
                                               1000 / 60 *
                                               exp(corr_preds_mamm)),
                               sdlog = my_sdlog)
    }
    x_mamm$prod <- x_mamm$p.ldd * x_mamm$mu
  } # End of "effective_dispersal" type


  # Otherwise, can just look at participation as seed dispersers
  if(type == "participation"){
    if(include_landscape_effect){
      x_bird$prod <- particip_preds_bird %*%
        t(plogis(-0.99 + 1.31 * coord_covars$treecover2000/100))
    } else{
      x_bird$prod <- particip_preds_bird
    }

    if(include_landscape_effect){
      x_mamm$prod <- particip_preds_mamm %*%
        t(plogis(-0.99 + 1.31 * coord_covars$treecover2000/100))
    } else{
      x_mamm$prod <- particip_preds_mamm
    }

  }

  # Now get these values for every site

  # First, only keep the species from x_bird and x_mamm present at each sites
  x_bird_by_site <- x_bird[which(x_bird$binomial %in% bird_presence_by_sitecoords_wide$binomial),]
  x_mamm_by_site <- x_mamm[which(x_mamm$binomial %in% mamm_presence_by_sitecoords_wide$binomial),]


  # Sum up the long-distance dispersal provided by each of these
  s_bird <- x_bird_by_site$prod * bird_presence_by_sitecoords_wide[,-1]
  s_mamm <- x_mamm_by_site$prod * mamm_presence_by_sitecoords_wide[,-1]


  # Counting how many species at each site
  dim(s_mamm)
  dim(s_bird)
  s_mamm_n <- s_mamm %>%
    summarise(across(1:ncol(s_mamm), function(x) sum(x > 0)))
  s_bird_n <- s_bird %>%
    summarise(across(1:ncol(s_bird), function(x) sum(x > 0)))
  s_both_n <- bind_rows(s_mamm, s_bird) %>%
    summarise(across(1:ncol(s_mamm), function(x) sum(x > 0)))

  # Version with total long-distance dispersal
  # First will get a version that is a sum across all species
  if(sum_or_mean == "sum"){
    disp_bird <- colSums(s_bird)
    disp_mamm <- colSums(s_mamm)
    disp_both <- disp_bird + disp_mamm
  } else{# Or can get a version that is the mean across the species
    disp_bird <- colSums(s_bird) / s_bird_n
    disp_mamm <- colSums(s_mamm) / s_mamm_n
    disp_both <- (colSums(s_bird) + colSums(s_mamm)) / s_both_n
  }






  # Output this
  output_list <- list(disp_bird, disp_mamm, disp_both)
  names(output_list) <- c("bird", "mamm", "both")
  return(output_list)

}
