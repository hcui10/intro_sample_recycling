## This file contains illustartion plots for nested Monte Carlo simulation accompanying 
## [Introducing Sample Recycling Method]
################################################################################

source("helpers_util.R", chdir = TRUE)
source("helpers_plot.R", chdir = TRUE)
source("algo.R", chdir = TRUE)

library(latex2exp)
library(ggplot2)

# Parameters for Illustrative Plot
F0 <- 100
N.out <- 2
N.in <- 3
tau <- 1/4
dt <- 1/252
t1 <- 4/52
mu <- 0.07
r <- 0.05 
d <- 0.02
sigma <- 0.3

# Helper Function: Simulate One Path
sim_path <- function(dt, tau, F0, mu, sigma) { 
  Wt <- cumsum(rnorm(as.integer(tau/dt), mean = 0, sd = sqrt(dt)))
  timeline <- seq(from = dt, to = tau, by = dt)
  return( c(F0, F0*exp( (mu-0.5*sigma^2)*timeline + sigma*Wt )) )
}

# Simulate Sample Paths
set.seed(234) # set seed for reproducibility 
Ft.samples <- # outer loop 
  sort(rlnorm(n = N.out, 
              meanlog = log(F0) + (mu - d - 0.5 * sigma^2) * t1, 
              sdlog   = sigma * sqrt(t1)))
Ft.sample.mat <- matrix(NA, 
                        nrow = N.out * N.in, 
                        ncol = as.integer( tau / dt ) + 1)
Ft.sample.mat <- # inner loop
  cbind(rep(F0, times = N.out * N.in), # append F0 upfront
        t(sapply(rep(Ft.samples, each = N.in), function(Ft) 
          sim_path(dt, tau, F0 = Ft, mu = r, sigma = sigma) )))
t.axis <- # time axis 
  c(0, t1, t1 + dt * seq_len(ncol(Ft.sample.mat)-2))
Ft.sample.df <- data.frame( # construct a data frame for plot
  "Time" = rep(t.axis, each = nrow(Ft.sample.mat)) * 52, 
  "Index" = rep(paste0("X", seq_len(nrow(Ft.sample.mat))), times = length(t.axis)), 
  "Value" = c(Ft.sample.mat), 
  "Group" = rep(rep(paste0("F", seq_len(N.out)), each = N.in), times = length(t.axis)), 
  "Inner_Flag" = c(rep("No", times = 2 * N.out * N.in), 
                   rep("Yes", times = (ncol(Ft.sample.mat) - 2) * N.out * N.in)))

# Plot the Illustrative Example
MC.concept.plot <- ggplot(data = Ft.sample.df, aes(x = Time)) + 
  # plot arrows
  geom_segment(aes(x = 0, xend = t1 * 52 * 0.97,
                   y = max(Ft.sample.mat), yend = max(Ft.sample.mat)),
               arrow = arrow(length = unit(0.20,"inches")), color = "grey") +
  geom_segment(aes(x = t1 * 52, xend = t.axis[length(t.axis)] * 52,
                   y = max(Ft.sample.mat), yend = max(Ft.sample.mat)),
               arrow = arrow(length = unit(0.20,"inches")), color = "grey") +
  geom_segment(aes(x = t1 * 52 * 0.97, xend = 0,
                   y = min(Ft.sample.mat) * 0.95, yend = min(Ft.sample.mat) * 0.95),
               arrow = arrow(length = unit(0.20,"inches")), color = "grey") +
  geom_segment(aes(x = t.axis[length(t.axis)] * 52, xend = t1 * 52,
                   y = min(Ft.sample.mat) * 0.95, yend = min(Ft.sample.mat) * 0.95),
               arrow = arrow(length = unit(0.20,"inches")), color = "grey") +
  # plot simulated sample paths
  geom_line(aes(y = Value, group = Index, color = Group, alpha = Inner_Flag)) + 
  # plot vertical dashed lines 
  geom_vline(xintercept = 0, 
             linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = t1 * 52, 
             linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = t.axis[length(t.axis)] * 52, 
             linetype = "dashed", color = "grey") + 
  # add annotations to arrows 
  annotate(geom = "text", x = mean(c(0, t1)) * 52, y = max(Ft.sample.mat), 
           label = "Real-World Measure", vjust = -0.5, color = "grey", size = 4) + 
  annotate(geom = "text", x = mean(c(t1, t.axis[length(t.axis)])) * 52, y = max(Ft.sample.mat), 
           label = "Risk-Neutral Measure", vjust = -0.5, color = "grey", size = 4) + 
  annotate(geom = "text", x = 0, y = min(Ft.sample.mat) * 0.95, 
           label = "Risk~Measure~rho", hjust = +0.27, vjust = -0.8, color = "grey", size = 4, parse = TRUE) + 
  annotate(geom = "text", x = t1 * 52, y = min(Ft.sample.mat)* 0.95 , 
           label = "Expected Loss", vjust = -1, color = "grey", size = 4) + 
  annotate(geom = "text", x = t.axis[length(t.axis)] * 52, y = min(Ft.sample.mat) * 0.95, 
           label = "Loss on Path", hjust = +0.7, vjust = -0.8, color = "grey", size = 4) + 
  # add labels 
  labs(x = TeX("Time $t$ : From 0 to $T$"), 
       y = TeX("$F_t$ : Risk Factor Level at Time $t$"), 
       title = "Illustration of Nested Monte Carlo Simulation") + 
  # add tick mark labels 
  scale_x_continuous(breaks = c(0, t1 * 52, t.axis[length(t.axis)] * 52), 
                     labels = c("Current Time", "Time of Interest", "Horizon")) + 
  # add legend 
  scale_color_manual(name = "Outer Scenario", 
                     labels = c("Outer Scenario #1", "Outer Scenario #2"), 
                     values = c("darkblue", "darkorange")) +
  scale_alpha_manual(name = "Inner Scenario", 
                     labels = c("Inner Paths for Scenario #1", 
                                "Inner Paths for Scenario #2"), 
                     values = c("No" = 0.4, "Yes" = 0.4)) +
  # adjust plot theme 
  guides(color = guide_legend(nrow = 2), 
         alpha = guide_legend(nrow = 2, 
                              override.aes = list(colour = c("darkblue", "darkorange"))))
MC.concept.plot <- apply_theme(MC.concept.plot)

print(MC.concept.plot)
ggsave(filename = "../img/MC_concept.png", plot = MC.concept.plot, 
       width = 20, height = 10, units = "in", dpi = 300)
