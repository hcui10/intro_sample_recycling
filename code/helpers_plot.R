## This file contains helper functions to generate plots accompanying 
## [Introducing Sample Recycling Method]
################################################################################

# theme helper
apply_theme <- function(plot.obj) {
  plot.obj <- plot.obj + theme_bw() + 
    theme(panel.border = element_blank(), 
          plot.title = element_text(size = 16, face = "bold"),
          legend.position = "bottom", 
          legend.title = element_text(size = 16), 
          legend.text = element_text(size = 16), 
          axis.line = element_line(colour = "black"), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16)) 
  return(plot.obj)
}

# visualize a matrix for estimated losses
make_loss_plot <- function(loss.mat, F.axis, L.axis, method.ylab, plotfile.path, 
                           return.plot.obj = FALSE) {
  # construct data frame for plot
  loss.est.df <- data.frame(
    cbind(F.axis, 
          t(apply(loss.mat, MARGIN = 2, # apply by column
                  FUN = function(arr) 
                    c(mean(arr, na.rm = TRUE), 
                      median(arr, na.rm = TRUE), 
                      sqrt(var(arr, na.rm = TRUE)), 
                      min(arr, na.rm = TRUE), 
                      max(arr, na.rm = TRUE))
          ))))
  rownames(loss.est.df) <- c()
  colnames(loss.est.df) <- 
    c("F_axis", "L_mean", "L_median", "L_sd", "L_min", "L_max")
  loss.est.df["True_Val"] <- L.axis
  # plot
  loss.est.plot <- ggplot(loss.est.df, aes(x = F_axis)) + 
    geom_ribbon(aes(ymin = L_min, ymax = L_max, 
                    fill = "violet"), alpha = 0.5) + 
    geom_ribbon(aes(ymin = L_mean - 3 * L_sd,
                    ymax = L_mean + 3 * L_sd, 
                    fill = "skyblue"), alpha = 0.5) +
    geom_line(aes(y = L_mean, color = "Mean_Est")) + 
    geom_line(aes(y = L_median, color = "Median_Est")) + 
    geom_line(aes(y = True_Val, color = "True_Val")) + 
    scale_fill_identity(name = "Volatility", guide = "legend", 
                        labels = c("Mean +/- 3 SD", "Min to Max")) + 
    scale_colour_manual(name = "Center of Estimates", 
                        values = c("True_Val" = "darkblue", 
                                   "Median_Est" = "greenyellow", 
                                   "Mean_Est" = "darkorange"), 
                        labels = c("Mean Estimate", 
                                   "Median Estimate", 
                                   "Theoretical Value")) + 
    labs(x = TeX("$F_t$ : Risk Factor Level at Time $t$"), 
         y = TeX("$\\widehat{L}_t(F_t)$"), 
         title = TeX(paste0("$\\widehat{L}_t(F_t)$", " Estimated Using ", method.ylab))) + 
    # coord_cartesian(
    #   ylim = c(min(loss.est.df[, c("True_Val", "L_mean", "L_median")]) * 0.9, 
    #            max(loss.est.df[, c("True_Val", "L_mean", "L_median")]) * 1.1)) + 
    guides(fill = guide_legend(nrow = 1), col = guide_legend(nrow = 1))
  loss.est.plot <- apply_theme(loss.est.plot)
  # save plot 
  if (!return.plot.obj) {
    print(loss.est.plot)
    ggsave(filename = plotfile.path, plot = loss.est.plot, 
           width = 18, height = 6, units = "in", dpi = 300)
  } else {
    return(list("plot" = loss.est.plot, "df" = loss.est.df))
  }
}

# visualize a comparison of multiple estimates
make_comp_plot <- function(loss.mats, F.axis, L.axis, plotfile.path, 
                           center.measures = c("Mean", "Median"), 
                           return.plot.obj = FALSE) {
  # construct matrices for plot
  loss.ests <- 
    lapply(loss.mats, function(loss.mat) 
      t(apply(loss.mat, MARGIN = 2, # apply by column 
              FUN = function(arr) c(mean(arr, na.rm = TRUE), median(arr, na.rm = TRUE))
      )))
  loss.ests.df <- data.frame(do.call("rbind", loss.ests))
  rownames(loss.ests.df) <- c()
  colnames(loss.ests.df) <- c("Mean", "Median")
  loss.ests.df[, "F_axis"] <- rep(F.axis, times = length(loss.mats))
  loss.ests.df[, "Method"] <- rep(names(loss.mats), each = length(F.axis))
  loss.true.df <- data.frame(cbind(F.axis, L.axis))
  rownames(loss.true.df) <- c()
  colnames(loss.true.df) <- c("F_axis", "True_Val")
  # plot
  loss.comp.plot <- ggplot(data = loss.ests.df, aes(x = F_axis))
  if ("Mean" %in% center.measures) {
    loss.comp.plot <- loss.comp.plot + 
      geom_line(aes(y = Mean, group = Method, colour = Method, linetype = "Mean"))
  }
  if ("Median" %in% center.measures) {
    loss.comp.plot <- loss.comp.plot + 
      geom_line(aes(y = Median, group = Method, colour = Method, linetype = "Median"))
  }
  loss.comp.plot <- loss.comp.plot + 
    geom_line(aes(y = True_Val, color = "Theoretical Value"), data = loss.true.df) + 
    labs(x = TeX("$F_t$ : Risk Factor Level at Time $t$"), 
         y = TeX("$\\widehat{L}_t(F_t)$ : Loss Estimates"), 
         linetype = "Center of Estimates") + 
    guides(linetype = guide_legend(nrow = 1), col = guide_legend(nrow = 1))
  loss.comp.plot <- apply_theme(loss.comp.plot)
  # save plot 
  if (!return.plot.obj) {
    print(loss.comp.plot)
    ggsave(filename = plotfile.path, plot = loss.comp.plot, 
           width = 18, height = 6, units = "in", dpi = 300)
  } else {
    return(list("plot" = loss.comp.plot, "df" = loss.ests.df))
  }
}

# visualize a comparison of errors of multiple estimates
make_comp_err_plot <- function(est.mats, x.axis, y.axis, plotfile.path, symbol.latex, 
                               error.measures = c("MAPE", "RMSE"), 
                               x.axis.percentage = FALSE, 
                               smooth.param = 0.1, 
                               return.plot.obj = FALSE) {
  # construct matrices for plot
  est.summary.list <- 
    lapply(est.mats, function(est.mat) 
      t(sapply(seq_len(ncol(est.mat)), 
               function(i) c("MAPE (%)" = MAPE(est.mat[,i], y.axis[i], na.rm = TRUE), 
                             "RMSE (%)" = RMSE(est.mat[,i], y.axis[i], na.rm = TRUE) / y.axis[i])
      )))
  est.summary.df <- data.frame(do.call("rbind", est.summary.list))
  rownames(est.summary.df) <- c()
  colnames(est.summary.df) <- c("MAPE", "RMSE")
  est.summary.df[, "x_axis"] <- rep(x.axis, times = length(est.mats))
  est.summary.df[, "Method"] <- rep(names(est.mats), each = length(x.axis))
  # plot
  comp.err.plot <- ggplot(data = est.summary.df, aes(x = x_axis))
  if (!is.null(smooth.param)) {
    if ("MAPE" %in% error.measures) {
      comp.err.plot <- comp.err.plot + 
        geom_smooth(aes(y = MAPE, group = Method, color = Method, fill = Method, linetype = "MAPE (%)"), 
                    method = "loess", formula = y ~ x, span = smooth.param)
    }
    if ("RMSE" %in% error.measures) {
      comp.err.plot <- comp.err.plot + 
        geom_smooth(aes(y = RMSE, group = Method, color = Method, fill = Method, linetype = "RMSE (%)"), 
                    method = "loess", formula = y ~ x, span = smooth.param)
    }
  } else {
    if ("MAPE" %in% error.measures) {
      comp.err.plot <- comp.err.plot + 
        geom_line(aes(y = MAPE, group = Method, color = Method, linetype = "MAPE (%)"))
    }
    if ("RMSE" %in% error.measures) {
      comp.err.plot <- comp.err.plot + 
        geom_line(aes(y = RMSE, group = Method, color = Method, linetype = "RMSE (%)"))
    }
  }
  if (x.axis.percentage) {
    comp.err.plot <- comp.err.plot + 
      scale_x_continuous(labels = scales::percent)
  }
  comp.err.plot <- comp.err.plot + 
    scale_y_continuous(labels = scales::percent) + 
    labs(x = TeX("$F_t$ : Risk Factor Level at Time $t$"), y = "Error (%)", 
         title = TeX(paste("Percentage (%) Error of", symbol.latex, sep = " ")), 
         linetype = "Error Measure") + 
    guides(linetype = guide_legend(nrow = 1, override.aes = list(fill = NA, color = "black")), 
           col = guide_legend(nrow = 1))
  comp.err.plot <- apply_theme(comp.err.plot)
  # save plot 
  if (!return.plot.obj) {
    print(comp.err.plot)
    ggsave(filename = plotfile.path, plot = comp.err.plot, 
           width = 18, height = 6, units = "in", dpi = 300)
  } else {
    return(list("plot" = comp.err.plot, "df" = est.summary.df))
  }
}

# Freedman–Diaconis rule to compute bin width / num bins for histogram
find_nbins_FDrule <- function(x.arr) {
  x.arr <- c(x.arr)
  x.arr <- x.arr[!is.na(x.arr)]
  quartiles <- as.numeric(quantile(x.arr, probs = c(0.25, 0.75)))
  IQR <- quartiles[2] - quartiles[1]
  n <- length(x.arr)
  # rule of thumb bin size 
  h <- 2 * IQR / (n^(1/3))
  nbins <- round( ( max(x.arr) - min(x.arr) ) / h , digits = 0)
  return(nbins)
}

# visualize histogram
make_hist <- function(df, nbins, true_val, return.plot.obj = FALSE, plotfile.path) {
  hist_plot <- ggplot(tidyr::gather(df, key = var_name, value = est_val), 
                      aes(x = est_val, fill = var_name, color = var_name)) + 
    geom_histogram(position = "identity", alpha = 0.5, bins = nbins) + 
    geom_vline(xintercept = true_val, linetype = "dashed") + 
    labs(x = "VaR at 95 % : Estimated Risk Measures", y = "Count", 
         fill = "Estimation Method", color = "Estimation Method", 
         title = TeX(paste0("Histogram of Estimated Risk Measures $\\widehat{\\rho} ( \\widehat{L}_t(F_t) )$ Based on Estimated Losses")))
  hist_plot <- apply_theme(hist_plot)
  # save plot 
  if (!return.plot.obj) {
    print(hist_plot)
    ggsave(filename = plotfile.path, plot = hist_plot, 
           width = 18, height = 6, units = "in", dpi = 300)
  } else {
    return(hist_plot)
  }
}
