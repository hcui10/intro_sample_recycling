## This file generates plots for 
## [Introducing Sample Recycling Method]
################################################################################

source("helpers_util.R", chdir = TRUE)
source("helpers_plot.R", chdir = TRUE)
source("algo.R", chdir = TRUE)

library(latex2exp)
library(ggplot2)

################################################################################

### INNER LAYER : Loss Est ### 
# Monte Carlo 
loss.est.MC <- readRDS("../data/inner_asian_result_MC.rds")[["loss_est"]]
sim.params <- readRDS("../data/inner_asian_result_MC_sim_params.rds")
loss_func <- get_loss_func(sim.params$inner.params, sim.params$contract.params)
make_loss_plot(loss.mat = loss.est.MC,
               F.axis = sim.params$Ft.mat[1,], 
               L.axis = loss_func(sim.params$Ft.mat[1,]), 
               method.ylab = method2fullname("MC"),
               plotfile.path = "../img/asian_loss_est_MC.png")
# Sample Recycling Method
loss.est.SRM <- readRDS("../data/inner_asian_result_recy_2.rds")[["loss_est"]]
make_loss_plot(loss.mat = loss.est.SRM,
               F.axis = sim.params$Ft.mat[1,], 
               L.axis = loss_func(sim.params$Ft.mat[1,]), 
               method.ylab = method2fullname("lambda_true"),
               plotfile.path = "../img/asian_loss_est_recy.png")
# Sample Recycling Method with Density Ratio Estimation 
loss.est.DRE <- readRDS("../data/inner_asian_result_recy_naive_2.rds")[["loss_est"]]
make_loss_plot(loss.mat = loss.est.DRE,
               F.axis = sim.params$Ft.mat[1,], 
               L.axis = loss_func(sim.params$Ft.mat[1,]), 
               method.ylab = method2fullname("lambda_naive"),
               plotfile.path = "../img/asian_loss_est_recy_naive.png")
# loss comparison 
Lt.est.ls <- list(loss.est.MC, loss.est.SRM, loss.est.DRE)
names(Lt.est.ls) <- c(method2fullname("MC"), 
                      method2fullname("lambda_true"), 
                      method2fullname("lambda_naive"))
make_comp_err_plot(est.mats = Lt.est.ls, # MAPE with LOESS smoothing 
                   x.axis = sim.params$Ft.mat[1,], y.axis = loss_func(sim.params$Ft.mat[1,]), 
                   symbol.latex = "Loss Estimate $\\widehat{L}_t(F_t)$", 
                   error.measures = c("MAPE"), plotfile.path = "../img/asian_loss_est_comp_err_MAPE.png")
make_comp_err_plot(est.mats = Lt.est.ls, # RMSE and MAPE without LOESS smoothing 
                   x.axis = sim.params$Ft.mat[1,], y.axis = loss_func(sim.params$Ft.mat[1,]), 
                   symbol.latex = "Loss Estimate $\\widehat{L}_t(F_t)$", 
                   smooth.param = NULL, plotfile.path = "../img/asian_loss_est_comp_err.png")
###

### OUTER LAYER : Risk Measure Est ###
# Histogram of VaR 95 % est
sim.params <- readRDS("../data/outer_asian_result_MC_sim_params.rds")
VaR_func <- get_VaR_func(sim.params$outer.params, sim.params$inner.params, sim.params$contract.params)
loss_func <- get_loss_func(sim.params$inner.params, sim.params$contract.params)
Lt.est.ls <- list(
  readRDS("../data/outer_asian_result_MC.rds")[["loss_est"]], 
  readRDS("../data/outer_asian_result_recy_100.rds")[["loss_est"]], 
  readRDS("../data/outer_asian_result_recy_naive_100.rds")[["loss_est"]]
)
names(Lt.est.ls) <- c(method2fullname("MC"), 
                      method2fullname("lambda_true"), 
                      method2fullname("lambda_naive"))
Lt.est.ls[["Theory"]] <- loss_func(sim.params$Ft.mat)
VaR.est.ls <- 
  lapply(Lt.est.ls, function(Lt.est.mat) 
    do.call("rbind", parallel::mclapply(
      data.frame(t(Lt.est.mat)), 
      function(Lt.est.arr) risk_measure_est(Lt.est.arr, risk.type = "VaR", est.params = NULL), 
      mc.cores = parallel::detectCores(), mc.allow.recursive = FALSE)))
VaR95.est.df <- data.frame(sapply(VaR.est.ls, function(VaR.est.mat) VaR.est.mat[,95]), row.names = NULL)
colnames(VaR95.est.df) <- names(VaR.est.ls)
make_hist(VaR95.est.df[,names(VaR95.est.df) != "Theory"], 
          nbins = max(sapply(VaR95.est.df[,names(VaR95.est.df) != "Theory"], find_nbins_FDrule)), 
          true_val = VaR_func(0.95), 
          plotfile.path = "../img/asian_VaR95_est_hist.png")
sapply(VaR95.est.df, function(est.arr) summarize(est.arr, true.val = VaR_func(0.95)))
# comparison of MAPE of VaR est - Sample Recycling Method 
probs.axis <- seq(from = 0.01, to = 0.99, by = 0.01)
Lt.est.ls <- list(
  "5 Reference Points" = 
    readRDS("../data/outer_asian_result_recy_5.rds")[["loss_est"]], 
  "20 Reference Points" = 
    readRDS("../data/outer_asian_result_recy_20.rds")[["loss_est"]], 
  "50 Reference Points" = 
    readRDS("../data/outer_asian_result_recy_50.rds")[["loss_est"]], 
  "100 Reference Points" = 
    readRDS("../data/outer_asian_result_recy_100.rds")[["loss_est"]]
)
VaR.est.ls <- 
  lapply(Lt.est.ls, function(Lt.est.mat) 
    do.call("rbind", parallel::mclapply(
      data.frame(t(Lt.est.mat)), 
      function(Lt.est.arr) risk_measure_est(Lt.est.arr, risk.type = "VaR", est.params = NULL), 
      mc.cores = parallel::detectCores(), mc.allow.recursive = FALSE)))
make_comp_err_plot(est.mats = VaR.est.ls, # MAPE with LOESS smoothing 
                   x.axis = probs.axis, y.axis = VaR_func(probs.axis), 
                   symbol.latex = "Risk Measure VaR $\\widehat{\\rho} ( v_0^t \\cdot \\widehat{L}_t(F_t) )$", 
                   error.measures = c("MAPE"), x.axis.percentage = TRUE, 
                   plotfile.path = "../img/asian_recy_VaR_est_err_MAPE.png")
make_comp_err_plot(est.mats = VaR.est.ls, # RMSE and MAPE without LOESS smoothing 
                   x.axis = probs.axis, y.axis = VaR_func(probs.axis), 
                   symbol.latex = "Risk Measure VaR $\\widehat{\\rho} ( v_0^t \\cdot \\widehat{L}_t(F_t) )$", 
                   smooth.param = NULL, x.axis.percentage = TRUE, 
                   plotfile.path = "../img/asian_recy_VaR_est_err.png")
# comparison of MAPE of VaR est - Sample Recycling Method with Density Ratio Estimation 
Lt.est.ls <- list(
  "5 Reference Points" = 
    readRDS("../data/outer_asian_result_recy_naive_5.rds")[["loss_est"]], 
  "20 Reference Points" = 
    readRDS("../data/outer_asian_result_recy_naive_20.rds")[["loss_est"]], 
  "50 Reference Points" = 
    readRDS("../data/outer_asian_result_recy_naive_50.rds")[["loss_est"]], 
  "100 Reference Points" = 
    readRDS("../data/outer_asian_result_recy_naive_100.rds")[["loss_est"]]
)
VaR.est.ls <- 
  lapply(Lt.est.ls, function(Lt.est.mat) 
    do.call("rbind", parallel::mclapply(
      data.frame(t(Lt.est.mat)), 
      function(Lt.est.arr) risk_measure_est(Lt.est.arr, risk.type = "VaR", est.params = NULL), 
      mc.cores = parallel::detectCores(), mc.allow.recursive = FALSE)))
make_comp_err_plot(est.mats = VaR.est.ls, # MAPE with LOESS smoothing 
                   x.axis = probs.axis, y.axis = VaR_func(probs.axis), 
                   symbol.latex = "Risk Measure VaR $\\widehat{\\rho} ( v_0^t \\cdot \\widehat{L}_t(F_t) )$", 
                   error.measures = c("MAPE"), x.axis.percentage = TRUE, 
                   plotfile.path = "../img/asian_recy_naive_VaR_est_err_MAPE.png")
make_comp_err_plot(est.mats = VaR.est.ls, # RMSE and MAPE without LOESS smoothing 
                   x.axis = probs.axis, y.axis = VaR_func(probs.axis), 
                   symbol.latex = "Risk Measure VaR $\\widehat{\\rho} ( v_0^t \\cdot \\widehat{L}_t(F_t) )$", 
                   smooth.param = NULL, x.axis.percentage = TRUE, 
                   plotfile.path = "../img/asian_recy_naive_VaR_est_err.png")
### 

### USE CASE ###
sim.params <- readRDS("../data/use_case_asian_result_recy_5_sim_params.rds")
VaR_func <- get_VaR_func(sim.params$outer.params, sim.params$inner.params, sim.params$contract.params)
loss_func <- get_loss_func(sim.params$inner.params, sim.params$contract.params)
set.seed(999)
idx.shuffle <- sample(seq_len(sim.params$outer.params$N.out), 
                      size = sim.params$outer.params$N.out, 
                      replace = FALSE)
Lt.est.ls <- list(
  c(readRDS("../data/use_case_asian_result_recy_5.rds")[["loss_est"]][1,idx.shuffle]), 
  c(readRDS("../data/use_case_asian_result_recy_naive_5.rds")[["loss_est"]][1,idx.shuffle])
)
names(Lt.est.ls) <- c(method2fullname("lambda_true"), method2fullname("lambda_naive"))
idx.ref <- which(sim.params$Ft.mat[1,idx.shuffle] %in% sim.params$Ft.ref.mat[1,])
idx.tar <- which(!(sim.params$Ft.mat[1,idx.shuffle] %in% sim.params$Ft.ref.mat[1,]))
MC.Lt.est <- Lt.est.ls[[method2fullname("lambda_true")]][idx.ref]
num.refs <- as.integer(seq(from = 50, to = 495, by = 5))
VaR95.est.df <- data.frame(sapply(Lt.est.ls, function(Lt.est.arr) sapply(num.refs, function(npts) as.numeric(
  quantile(c(MC.Lt.est, Lt.est.arr[idx.tar][1:npts]), probs = 0.95)))))
VaR95.est.df[["npts"]] <- num.refs
VaR95.est.plot.df <- tidyr::gather(VaR95.est.df, key = `Estimation Method`, value = est_val, -npts)
use_case_plot <- ggplot(VaR95.est.plot.df, aes(x = npts)) +
  geom_smooth(aes(y = est_val, group = `Estimation Method`, 
                  color = `Estimation Method`, fill = `Estimation Method`), 
              method = "loess", formula = y ~ x, span = 0.1) +
  geom_point(aes(y = est_val, group = `Estimation Method`, color = `Estimation Method`), 
             shape = 20, alpha = 0.5) + 
  geom_hline(yintercept = VaR_func(0.95), linetype = "dashed", color = "steelblue4") + 
  labs(x = TeX("Number of Additional $F_t$ Samples"), y = TeX("Estimated Risk Measure $\\widehat{\\rho}$"),
       title = "Risk Measure Estimation using Sample Recycling Methods") 
use_case_plot <- apply_theme(use_case_plot)
print(use_case_plot)
ggsave(filename = "../img/use_case.png", plot = use_case_plot, 
       width = 18, height = 6, units = "in", dpi = 300)
# print summary 
round(matrix(c(sort(sim.params$Ft.ref.mat), 
               loss_func(sort(sim.params$Ft.ref.mat)), 
               MC.Lt.est), 
             nrow = 3, byrow = TRUE), 
      digits = 2)
sapply(Lt.est.ls, function(Lt.est.arr) as.numeric(quantile(Lt.est.arr, probs = 0.95)))
###
