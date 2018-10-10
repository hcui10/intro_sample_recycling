## This file contains illustartion plots for density ratio estimation accompanying 
## [Introducing Sample Recycling Method]
################################################################################

source("helpers_util.R", chdir = TRUE)
source("algo.R", chdir = TRUE)

library(latex2exp)
library(ggplot2)

# set seed for reproducibility 
set.seed(999)

# simulate data
lambda_params <- GBM_lambda(F.tar = 110, F.ref = 105, sigma = 0.2, r = 0, d = 0, dt = 1)
x_tar <- rlnorm(n = 500, meanlog = lambda_params$tar_meanlog, sdlog = lambda_params$tar_sdlog)
x_ref <- rlnorm(n = 500, meanlog = lambda_params$ref_meanlog, sdlog = lambda_params$ref_sdlog)

# plot histogram
hist.breaks <- # construct break points 
  c(min(x_tar, x_ref), 
    as.numeric(head( # remove first and last (sample min and max) from quantiles 
      quantile(x_ref, probs = seq(from = 0, to = 1, length.out = 10 + 1)[-1]), -1)), 
    max(x_tar, x_ref))
hist.plot <- ggplot(data = data.frame(x_ref), aes(x = x_ref)) + 
  # generate histograms
  geom_histogram(aes(x = x_ref, y = ..density.., fill = "Reference"), breaks = hist.breaks, 
                 alpha = 0.3, position = "identity") + 
  geom_histogram(aes(x = x_tar, y = ..density.., fill = "Target"), breaks = hist.breaks, 
                 alpha = 0.3, position = "identity") + 
  # show counts 
  stat_bin(aes(x = x_ref, y = ..density.., label = ..count..), 
           geom = "text", vjust = -0.0, breaks = hist.breaks, color = "darkblue") + 
  stat_bin(aes(x = x_tar, y = ..density.., label = ..count..), 
           geom = "text", vjust = +0.0, breaks = hist.breaks, color = "darkorange") + 
  # overlay density functions 
  stat_function(fun = dlnorm, color = "darkorange", 
                args = list(meanlog = lambda_params$tar_meanlog, 
                            sdlog   = lambda_params$tar_sdlog)) + 
  stat_function(fun = dlnorm, color = "darkblue", 
                args = list(meanlog = lambda_params$ref_meanlog, 
                            sdlog   = lambda_params$ref_sdlog)) + 
  # visualize break points 
  geom_vline(xintercept = hist.breaks, linetype = "dashed", alpha = 0.2) + 
  # names axis
  labs(x = TeX("$F$ : Risk Factor Level"), y = "Density") + 
  # add legend
  scale_fill_manual(name = "Sample Distribution Given Starting Value", 
                    values = c("Reference" = "darkblue", "Target" = "darkorange"), 
                    labels = c("Reference", "Target")) + 
  # adjust theme and delete background and grids
  theme_bw() + theme(axis.line = element_line(colour = "black"), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.border = element_blank(), 
                     legend.position = "bottom") 
print(hist.plot)
ggsave(filename = "../img/ratio_est_hist.png", plot = hist.plot, 
       width = 7.5, height = 7.5, units = "in", dpi = 300)

# estimate density ratio with varying number of blocks
DR_est_10 <- lambda_step_approx(x_tar, x_ref, n_block = 10, x_min = 0, x_max = Inf)
# DR_est_15 <- lambda_step_approx(x_tar, x_ref, n_block = 15, x_min = 0, x_max = Inf)
# DR_est_20 <- lambda_step_approx(x_tar, x_ref, n_block = 20, x_min = 0, x_max = Inf)

# construct data frames for plot
x.axis <- seq(from = min(hist.breaks), to = max(hist.breaks), length.out = 100)
ratio.df <- data.frame(
  "x" = rep(x.axis, times = 2), 
  "ratio_pred" = c(lambda_params$lambda(x.axis), 
                   DR_est_10(x.axis)), 
                   # DR_est_15(x.axis), 
                   # DR_est_20(x.axis)), 
  "est_method" = c(rep("Theoretical Value"        , times = length(x.axis)), 
                   rep("Estimation with 10 Blocks", times = length(x.axis)))#, 
                   # rep("Estimation with 15 Blocks", times = length(x.axis)), 
                   # rep("Estimation with 20 Blocks", times = length(x.axis)))
)

# construct annotation labels
x.breaks <- which(diff(ratio.df[ratio.df$est_method == "Estimation with 10 Blocks","ratio_pred"]) != 0)
x.start  <- c(1, x.breaks + 1)
x.end    <- c(x.breaks, length(x.axis))
x.label  <- vector(mode = "numeric", length = 10)
for (i in seq_len(10)) x.label[i] <- median(x.axis[x.start[i]:x.end[i]])
y.label <- DR_est_10(x.label)

# plot density ratio estimation
ratio.plot <- ggplot(data = ratio.df, aes(x = x, y = ratio_pred, color = est_method)) + 
  # geom_step(data = subset(ratio.df, est_method != "Theoretical Value")) + 
  geom_line(data = subset(ratio.df, est_method == "Theoretical Value"), lwd = 1) + 
  geom_step(data = subset(ratio.df, est_method == "Estimation with 10 Blocks"), lwd = 1) + 
  annotate(geom = "text", x = x.label, y = y.label, label = sprintf("%.2f", round(y.label, digits = 2)), 
           hjust = 0.4, vjust = -0.3, size = 4) + 
  geom_vline(xintercept = hist.breaks, linetype = "dashed", alpha = 0.2) + 
  labs(x = TeX("$F$ : Risk Factor Level"), 
       y = TeX("$\\widehat{\\lambda}(F)$ : Estimated Density Ratios"), 
       color = "Number of Blocks") + 
  # guides(col = guide_legend(nrow = 2)) + 
  # adjust theme and delete background and grids
  theme_bw() + theme(axis.line = element_line(colour = "black"), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.border = element_blank(), 
                     legend.position = "bottom") 
print(ratio.plot)
ggsave(filename = "../img/ratio_est_ratio.png", plot = ratio.plot, 
       width = 7.5, height = 7.5, units = "in", dpi = 300)
