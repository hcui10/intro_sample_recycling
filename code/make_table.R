## This file contains helper functions to generate tables accompanying 
## [Introducing Sample Recycling Method]
################################################################################

source("helpers_util.R", chdir = TRUE)
source("helpers_tab.R", chdir = TRUE)
source("algo.R", chdir = TRUE)

################################################################################

sim.params <- readRDS("../data/inner_asian_result_MC_sim_params.rds")
sim.results <- list(
  "MC"           = readRDS("../data/inner_asian_result_MC.rds"), 
  "lambda_true"  = readRDS("../data/inner_asian_result_recy_2.rds"), 
  "lambda_naive" = readRDS("../data/inner_asian_result_recy_naive_2.rds")
)
loss_func <- get_loss_func(sim.params$inner.params, sim.params$contract.params)

# data validation 
cat("\n")
check_data(sim.results, F.axis = sim.params$Ft.mat[1,], F.ref.axis = sim.params$Ft.ref.mat[1,])
cat("\n")

# print parameters
print(paste(paste(rep("#", times = 20), collapse = ""), "CONTRACT PARAMETERS", 
            paste(rep("#", times = 20), collapse = ""), sep = " "))
print(t(data.frame(lapply(
  sim.params$contract.params, FUN = function(x) ifelse(is.null(x), "NULL", as.character(x))))))
print(paste(paste(rep("#", times = 20), collapse = ""), "SIMULATION PARAMETERS", 
            paste(rep("#", times = 20), collapse = ""), sep = " "))
print(round(t(data.frame(sim.params$inner.params)), digits = 2))

# details of each estimation method
for (method.type in sim.params$sim.methods) {
  # get method full name
  method.full.name <- method2fullname(method.type)
  print(paste(paste(rep("-", times = 20), collapse = ""), 
              method.full.name, 
              paste(rep("-", times = 20), collapse = ""), sep = " "))
  cat("\n")
  # time profile
  print("TIME PROFILE in Milli-Seconds")
  time.profile.mat <- make_time_table(method.type, sim.results, 
                                      F.axis = sim.params$Ft.mat[1,],
                                      F.ref.axis = sim.params$Ft.ref.mat[1,])
  print(round(1000 * time.profile.mat, digits = 2))
  cat("\n")
}
