## This file contains illustartion plots for sample recycling method accompanying 
## [Introducing Sample Recycling Method]
################################################################################

source("helpers_util.R", chdir = TRUE)
source("helpers_plot.R", chdir = TRUE)
source("algo.R", chdir = TRUE)

library(latex2exp)
library(ggplot2)

F.ref <- 100
F.tar <- 98
dt <- 1/52
tau <-dt * 10

set.seed(123)
F.ref.mat <- matrix(NA, nrow = 3, ncol = as.integer(tau/dt) + 1)
F.ref.mat[1,] <- sim_path(dt, tau, F0 = F.ref, mu = 0.05, sigma = 0.3)
F.ref.mat[2,] <- sim_path(dt, tau, F0 = F.ref, mu = 0.05, sigma = 0.3)
F.ref.mat[3,] <- sim_path(dt, tau, F0 = F.ref, mu = 0.05, sigma = 0.3)
F.ref.mat <- t(F.ref.mat)

F.tar.mat <- F.ref.mat * 0.997
F.tar.mat[1,] <- rep(F.tar, times = ncol(F.ref.mat))
Ft.mat <- t(cbind(F.ref.mat, F.tar.mat))
t.axis <- seq_len(ncol(Ft.mat)) - 1
t.axis <- t.axis + 2
t.axis[1] <- 0
Ft.df <- data.frame(
  "Time" = rep(t.axis, each = nrow(Ft.mat)), 
  "Index" = rep(paste0("X", seq_len(nrow(Ft.mat))), times = length(t.axis)), 
  "Value" = c(Ft.mat), 
  "Group" = rep(rep(paste0("F", seq_len(2)), each = ncol(F.ref.mat)), times = length(t.axis)))

SR.concept.plot <- ggplot() + 
  geom_segment(aes(x = t.axis[1], xend = t.axis[length(t.axis)],
                   y = max(Ft.mat) * 1.005, yend = max(Ft.mat) * 1.005),
               arrow = arrow(length = unit(0.20,"inches")), color = "grey") + 
  annotate(geom = "text", x = mean(c(t.axis[1], t.axis[length(t.axis)])), y = max(Ft.mat) * 1.005, 
           label = "Risk-Neutral Measure", vjust = -0.5, color = "grey", size = 5) + 
  geom_vline(xintercept = t.axis[1], 
             linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = t.axis[2], 
             linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = t.axis[length(t.axis)], 
             linetype = "dashed", color = "grey") + 
  geom_line(aes(x = Time, y = Value, group = Index, linetype = Group, color = Group), data = Ft.df) + 
  geom_point(aes(x = rep(t.axis[2], times = ncol(F.ref.mat)), 
                 y = c(F.ref.mat[2,])), size = 3, color = "darkblue") + 
  geom_point(aes(x = 0, y = unique(F.ref.mat[1,])), size = 3, color = "darkblue") + 
  geom_point(aes(x = 0, y = unique(F.tar.mat[1,])), size = 3, color = "darkorange") + 
  coord_cartesian(ylim = c(min(Ft.mat) * 0.99, max(Ft.mat) * 1.01)) + 
  labs(x = TeX("Time : From $t$ to $T$"), 
       y = TeX("$F_t$ : Risk Factor Level at Time $t$"), 
       title = "Illustration of Sample Recycling Method for Inner Loop Simulation") + 
  scale_color_manual(name = "Outer Scenario", 
                     labels = c("Reference Scenario", "Target Scenario"), 
                     values = c("darkblue", "darkorange")) +
  scale_linetype_manual(name = "Outer Scenario", 
                        labels = c("Reference Scenario", "Target Scenario"), 
                        values = c("solid", "dashed")) + 
  scale_x_continuous(breaks = c(t.axis[1], t.axis[2], t.axis[length(t.axis)]), 
                     labels = c("Time of Interest : t", expression(t~+~Delta~t), "Horizon : T"))
SR.concept.plot <- apply_theme(SR.concept.plot)

print(SR.concept.plot)
ggsave(filename = "../img/SR_concept.png", plot = SR.concept.plot, 
       width = 20, height = 10, units = "in", dpi = 300)
