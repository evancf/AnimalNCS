# You need to run 04_displacement_bayes_model.R to the end with the terr_mod stuff
# as well as the 20_fragmentation_effect.R code

pres_dat <- sjPlot::get_model_data(tc_mod, type = "pred",
                                   terms="treecover2000 [all]")


dist_dat <- sjPlot::get_model_data(terr_mod_mean, type = "pred",
                                   terms = c("mean_fp [all]"))



# First presence

png(file = "./outputs/figures/both habitat effects.png",
    width = 6.5,
    height = 4.75,
    units = "in",
    res = 440)

par(mfrow=c(1,2),
    mar = c(5.1,4.1,2.1,2.1))
pres_xx <- 1-pres_dat$x
plot(NA,#pres_xx,
     NA,#pres_dat$predicted,
     xlab = "Landscape fragmentation",
     ylab = "Presence probability",
     xlim = c(0, 1),
     ylim = c(0,1),#c(0.2, 0.7),
     type = "l",
     frame = F,
     lwd = 2,
     las = 1,
     xaxt = "n")
axis(1, at = c(0,0.5,1))



site_summary <- timm_gfc %>%
  group_by(coord_id) %>%
  dplyr::summarize(avg = mean(Presabs),
                   frag = first(1-treecover2000))

points(avg ~ frag,
       data = site_summary,
       pch = 16,
       cex = 0.5,
       col = "grey"#,  rgb(0,0,0,0.05)
)


points(pres_xx,
       pres_dat$predicted,
       lwd = 2,
       type = "l")

text(-0.1,
     1,#0.7,
     "a", pos = 3, xpd = T,
     cex = 1.25, font = 2)

polygon(x = c(pres_xx, rev(pres_xx)),
        y = c(pres_dat$conf.low, rev(pres_dat$conf.high)),
        col = rgb(0,0,0,0.3),
        lty = 0)



# Next distance
dist_xx <- dist_dat$x
plot(NA,#dist_xx,
     NA,#dist_dat$predicted,
     xlab = "Human footprint index",
     ylab = "Animal displacement (m)",
     xlim = c(0, 1),
     ylim = c(0, 3000),
     type = "l",
     frame = F,
     lwd = 2,
     las = 1,
     xaxt = "n")
axis(1, at = c(0,0.5,1),
     labels = c(0, 25, 50)) # Back to the original scale




set.seed(1)
plot_df <- terr_df2[sample(1:nrow(terr_df2)),] %>%
  group_by(frug_event_id) %>%
  slice(1:3)

points(ifelse(plot_df$displacement > 3000, NA, plot_df$displacement)  ~ plot_df$mean_fp,
       data = plot_df,
       ylim = c(1, 3000),
       pch = 16,
       cex = 0.2,
       col = "grey"#,  rgb(0,0,0,0.05)
)

hist_df <- plot_df %>%
  filter(displacement > 3000)

dens_vals <- density(hist_df$mean_fp) %>%
  as.data.frame() %>%
  filter(x > 0, x < 1)

polygon(x = c(1,0,dens_vals$x),
        y = c(0,0,dens_vals$y) * 100 + 3050,
        col = "grey",
        border = NA,
        xpd = T)




points(dist_xx,
       dist_dat$predicted,
       type = "l",
       lwd = 2)

text(-0.1,3000, "b", pos = 3, xpd = T,
     cex = 1.25, font = 2)

polygon(x = c(dist_xx, rev(dist_xx)),
        y = c(dist_dat$conf.low, rev(dist_dat$conf.high)),
        col = rgb(0,0,0,0.3),
        lty = 0,
        xpd = T)




dev.off()

