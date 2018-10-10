## This file contains code for experiments covered in  
## [Introducing Sample Recycling Method]
################################################################################

# source scripts 
source("helpers_util.R", chdir = TRUE)
source("algo.R", chdir = TRUE)
# removes all objects except for functions
rm(list = setdiff(ls(), lsf.str()))

################################################################################
###### Helper Function for Inner Layer Estimate: Wrap Functions in algo.R ######
################################################################################
run_inner_layer <- function(contract.params, inner.params, 
                            method.type, Ft.mat, Ft.ref.mat, 
                            save.flag = TRUE, 
                            save.to.file = "../data/temp/L_est.rds", 
                            ncpus = 1, verbose = TRUE) {
  # unpacking params
  func.envir <- list2env(contract.params, envir = environment())
  func.envir <- list2env(inner.params, envir = environment())
  # initialize a storage list for simulation results
  sim.results.method <- list()
  if (method.type == "MC") {
    # Monte Carlo Simulation
    sim.results.method[["loss_est"]] <-
      matrix(NA, nrow = nrow(Ft.mat), ncol = ncol(Ft.mat))
    sim.results.method[["time"]][["sim"]] <-
      sim.results.method[["time"]][["eval"]] <- 
      matrix(0, nrow = nrow(Ft.mat), ncol = ncol(Ft.mat))
    for (i in seq_len(nrow(Ft.mat))) {
      for (j in seq_len(ncol(Ft.mat))) {
        est.MC <-
          option_price_MC(S0 = Ft.mat[i,j], K = K, sigma = sigma, r = r, d = d, tau = tau,
                          loss.type = loss.type, option.type = option.type,
                          num.MC = N.in, dt = dt,
                          # redundant parameters for European options
                          avg.target = avg.target, avg.method = avg.method, avg.step = 1,
                          ncpus = ncpus, timeit = TRUE)
        # work with L'Ecuyer-CMRG RNG to ensure reproducibility 
        skip_RNG_streams(ncpus, envir = globalenv()) 
        # record simulation results
        sim.results.method[["loss_est"]][i,j] <- est.MC$est
        sim.results.method[["time"]][["sim"]][i,j] <- 
          est.MC$sim.time["elapsed"]
        sim.results.method[["time"]][["eval"]][i,j] <- 
          est.MC$loss.eval.time["elapsed"]
      }
      # print progress
      if (verbose) {
        if ( i == 1 ) {
          cat("\n")
          print(paste0("MONTE CARLO PROGRESS: ", i, " / ", N.sim, " Simulations", 
                       " AT ", as.character(Sys.time())))
          print(sim.results.method[["loss_est"]][1,])
          cat("\n")
        } else if ( i %% 5 == 0 ) {
          print(paste0("MONTE CARLO PROGRESS: ", i, " / ", N.sim, " Simulations", 
                       " AT ", as.character(Sys.time())))
        }
      }
    }
  } else {
    classifier.type <- method2classifier(method.type)
    # Sample Recycling Method (potentially with Density Ratio Estimation)
    sim.results.method[["loss_est"]] <- 
      sim.results.method[["exp_lambda"]] <- 
      matrix(NA, nrow = nrow(Ft.mat), ncol = ncol(Ft.mat))
    sim.results.method[["time"]][["sim"]] <-
      sim.results.method[["time"]][["eval"]] <-
      sim.results.method[["time"]][["ratio_sim"]] <-
      sim.results.method[["time"]][["ratio_est"]] <-
      sim.results.method[["time"]][["ratio_eval"]] <-
      matrix(0, nrow = nrow(Ft.mat), ncol = ncol(Ft.mat))
    for (i in seq_len(nrow(Ft.mat))) {
      # reference - target mapping
      idx.tars <- find_tar_idx(Ft.ref.mat[i,], Ft.mat[i,])
      for (F.ref in Ft.ref.mat[i,]) {
        # get index of reference point in F.axis
        idx.ref <- which(Ft.mat[i,] == F.ref) 
        # Monte Carlo for reference point
        est.ref <-
          option_price_MC(S0 = F.ref, K = K, sigma = sigma, r = r, d = d, tau = tau,
                          loss.type = loss.type, option.type = option.type,
                          num.MC = N.ref, dt = dt,
                          # redundant parameters for European options
                          avg.target = avg.target, avg.method = avg.method, avg.step = 1,
                          ncpus = ncpus, timeit = TRUE)
        # work with L'Ecuyer-CMRG RNG to ensure reproducibility
        skip_RNG_streams(ncpus, envir = globalenv())
        # record simulation results 
        sim.results.method[["loss_est"]][i,idx.ref] <- est.ref$est
        sim.results.method[["time"]][["sim"]][i,idx.ref] <-
          est.ref$sim.time["elapsed"]
        sim.results.method[["time"]][["eval"]][i,idx.ref] <-
          est.ref$loss.eval.time["elapsed"]
        # extract reference point info for target points' use
        ref.paths.next.step <- est.ref$samples[2,]
        ref.losses <- est.ref$disc.losses
        # simulate reference sample for density ratio estimation
        if (method.type != "lambda_true") {
          sim.results.method[["time"]][["ratio_sim"]][i,idx.ref] <-
            system.time({
              x_de <- rlnorm(n = N.est, 
                             meanlog = log(F.ref) + (r - d - 0.5 * sigma^2) * dt, 
                             sdlog = sigma * sqrt(dt))
            })["elapsed"]
        } else { x_de <- NULL }
        # for each target point under current reference point
        idx.tar <- idx.tars[[toString(idx.ref)]]
        if ( any(is.na(idx.tar)) ) next 
        tar.results <- simplify2array( # each col is a one iteration
          parallel::mclapply(idx.tar, FUN = function(j) {
            if (method.type != "lambda_true") {
              # simulate target sample for density ratio estimation
              time.ratio.sim <- as.numeric(system.time({
                x_nu <- rlnorm(n = N.est,
                               meanlog = log(Ft.mat[i,j]) + (r - d - 0.5 * sigma^2) * dt,
                               sdlog = sigma * sqrt(dt))
              })["elapsed"])
            } else {
              x_nu <- NULL
              time.ratio.sim <- 0
            }
            # fit density ratios
            classifier.params <- switch(
              method.type,
              "lambda_true" = list("F.tar" = Ft.mat[i,j], "F.ref" = F.ref,
                                   "sigma" = sigma, "r" = r, "d" = d, "dt" = dt),
              "lambda_naive" = list("n_block" = 10, "x_min" = 0, "x_max" = Inf))
            time.ratio.est <- as.numeric(system.time({
              DR.ratio.est <-
                ratio_est(classifier.type = classifier.type, x_nu, x_de,
                          params = classifier.params)
            })["elapsed"])
            # estimate density ratios
            time.ratio.eval <- as.numeric(system.time({
              ratio.pred <- as.numeric(DR.ratio.est(ref.paths.next.step))
            })["elapsed"])
            # estimated density ratio quality
            est.exp.lambda <-
              mean(ratio.pred, na.rm = TRUE)
            # estimate loss
            time.eval <- as.numeric(system.time({
              est.loss <-
                mean(ratio.pred * ref.losses, na.rm = TRUE)
            })["elapsed"])
            # return
            return(c("time.ratio.sim" = time.ratio.sim,
                     "time.ratio.est" = time.ratio.est,
                     "time.ratio.eval" = time.ratio.eval,
                     "est.exp.lambda" = est.exp.lambda,
                     "time.eval" = time.eval,
                     "est.loss" = est.loss))
          }, mc.cores = ncpus, mc.allow.recursive = FALSE))
        # work with L'Ecuyer-CMRG RNG to ensure reproducibility
        skip_RNG_streams(ncpus, envir = globalenv())
        # assign values to storage matrices
        sim.results.method[["loss_est"]][i, idx.tar] <-
          tar.results["est.loss", ]
        sim.results.method[["exp_lambda"]][i, idx.tar] <-
          tar.results["est.exp.lambda", ]
        sim.results.method[["time"]][["eval"]][i, idx.tar] <-
          tar.results["time.eval", ]
        sim.results.method[["time"]][["ratio_sim"]][i, idx.tar] <-
          tar.results["time.ratio.sim", ]
        sim.results.method[["time"]][["ratio_est"]][i, idx.tar] <-
          tar.results["time.ratio.est", ]
        sim.results.method[["time"]][["ratio_eval"]][i, idx.tar] <-
          tar.results["time.ratio.eval", ]
      }
      # print progress
      if (verbose) {
        if ( i == 1 ) {
          cat("\n")
          print(paste0("SAMPLE RECYCLING WITH ", toupper(method.type), 
                       " PROGRESS: ", i, " / ", N.sim, " Simulations", 
                       " AT ", as.character(Sys.time())))
          print(sim.results.method[["loss_est"]][1,])
          cat("\n")
        } else if( i %% 10 == 0) {
          print(paste0("SAMPLE RECYCLING WITH ", toupper(method.type), 
                       " PROGRESS: ", i, " / ", N.sim, " Simulations", 
                       " AT ", as.character(Sys.time())))
        }
      }
    }
  }
  # save / return estimation results 
  if (save.flag) saveRDS(sim.results.method, file = save.to.file)
  else return(sim.results.method)
}
################################################################################

# ################################################################################
# ############################ European Option Params ############################
# ################################################################################
# contract.params.euro <- 
#   list("tau" = 0.25, "K" = 100, "sigma" = 0.3, "r" = 0.05, "d" = 0.02, 
#        "loss.type"  = "European", "option.type" = "call", 
#        "avg.target" = NULL, "avg.method" = NULL) # redundant params for European option
# inner.params.euro <- 
#   list("N.sim" = 100 ,   # run 100 simulations for each Ft value
#        "N.in"  = 1000,   # number of inner paths per scenario for Monte Carlo
#        "N.ref" = 5000,   # number of reference sample
#        "N.est" = 500 ,   # number of samples for density ratio estimation
#        "dt"    = 0.25    # tau / n.dt
#   )
# inner.params.euro[["n.dt"]] <- # tau / dt = 1, path independent
#   as.integer(contract.params.euro$tau / inner.params.euro$dt)
# outer.params.euro <- 
#   list("F0" = 100, "t1" = 1/253, "mu" = 0.07, "N.out" = 500, "num.ref" = 5)
# ################################################################################

################################################################################
############################# Asian Option Params ##############################
################################################################################
contract.params.asian <- 
  list("tau" = 0.25, "K" = 100, "sigma" = 0.3, "r" = 0.05, "d" = 0.02, 
       "loss.type"  = "Asian", "option.type" = "call", 
       "avg.target" = "price", "avg.method"  = "geometric")
inner.params.asian <- 
  list("N.sim" = 100 ,   # run 100 simulations for each Ft value
       "N.in"  = 1000,   # number of inner paths per scenario for Monte Carlo
       "N.ref" = 5000,   # number of reference sample
       "N.est" = 500 ,   # number of samples for density ratio estimation
       "dt"    = 1/52    # tau / n.dt
       )
inner.params.asian[["n.dt"]] <- # tau / dt = 13, path dependent
  as.integer(contract.params.asian$tau / inner.params.asian$dt)
outer.params.asian <- 
  list("F0" = 100, "t1" = 1/253, "mu" = 0.07, 
       "N.out" = 500, # consider 500 outer scenarios
       "num.ref" = 5) # with 5 reference scenarios for sample recycling methods
################################################################################

# print number of cores to be used
ncpus <- parallel::detectCores()
print(paste0("NUMBER OF CORES TO BE USED: ", ncpus))

# reproducible random numbers in parallel paradigm
RNGkind("L'Ecuyer-CMRG") 

################################################################################
##################### Estimate Loss at Time t with Fixed Ft ####################
################################################################################
##### Asian Option #####
# set up parameters
sim.params <- list(
  "contract.params" = contract.params.asian,
  "inner.params"    = inner.params.asian,
  "outer.params"    = list("N.out" = 500, "num.ref" = 2),
  "sim.methods"     = c("MC", "lambda_true", "lambda_naive")
)
sim.params[["Ft.mat"]] <- matrix(
  rep(seq(from = 90, to = 110,
          length.out = sim.params$outer.params$N.out),  # N.out columns
      times = sim.params$inner.params$N.sim),           # N.sim rows
  nrow = sim.params$inner.params$N.sim, byrow = TRUE)
sim.params[["Ft.ref.mat"]] <-
  t(apply(sim.params[["Ft.mat"]], MARGIN = 1, # apply by row (by each replication)
          FUN = function(row) as.numeric(
            quantile(row, probs =
                       seq_len(sim.params$outer.params$num.ref) /
                       ( sim.params$outer.params$num.ref + 1 ),
                     type = 1))))
# run experiments with different methods
set.seed(142857) # set seed for reproducibility
for ( method.type in sim.params$sim.methods ) {
  # report start
  cat("\n")
  print(paste0(paste(rep("#", times = 20), collapse = ""), " ",
               "RUNNING [FIXED Ft - ASIAN OPTION] WITH ",
               toupper(method2fullname(method.type)), " ",
               paste(rep("#", times = 20), collapse = "")))
  cat("\n")
  print(paste("START TO RUN", toupper(method2fullname(method.type)), "AT",
              as.character(Sys.time()), sep = " "))
  cat("\n")
  # run inner layer simulation
  run_inner_layer(
    sim.params$contract.params, sim.params$inner.params, method.type,
    sim.params$Ft.mat, sim.params$Ft.ref.mat,
    save.to.file = ifelse(method.type == "MC",
                          "../data/inner_asian_result_MC.rds",
                          paste0("../data/inner_asian_result_",
                                 method2filesuffix(method.type), "_",
                                 sim.params$outer.params$num.ref, ".rds")),
    ncpus = ncpus)
  saveRDS(sim.params,
          file = ifelse(method.type == "MC",
                        "../data/inner_asian_result_MC_sim_params.rds",
                        paste0("../data/inner_asian_result_",
                               method2filesuffix(method.type), "_",
                               sim.params$outer.params$num.ref,
                               "_sim_params.rds")))
  # report end
  print(paste("END RUNNING", toupper(method2fullname(method.type)), "AT",
              as.character(Sys.time()), sep = " "))
  cat("\n")
}
################################################################################

################################################################################
#################### Estimate Loss at Time t with Random Ft ####################
################################################################################
##### Asian Option #####
# set up parameters
sim.params <- list(
  "contract.params" = contract.params.asian,
  "inner.params"    = inner.params.asian,
  "outer.params"    = outer.params.asian,
  "sim.methods"     = c("MC", "lambda_true", "lambda_naive")
)
set.seed(285714) # set seed for reproducibility
Ft.RNG <- get_Ft_RNG(sim.params$outer.params, sim.params$contract.params)
sim.params[["Ft.mat"]] <- t(simplify2array(parallel::mclapply(
  seq_len(sim.params$inner.params$N.sim), FUN = function(i)
    sort(Ft.RNG(n = sim.params$outer.params$N.out)),
  mc.cores = ncpus, mc.allow.recursive = FALSE)))
skip_RNG_streams(ncpus, # work with L'Ecuyer-CMRG RNG to ensure reproducibility
                 envir = globalenv())
# run experiments with different methods
for ( method.type in sim.params$sim.methods ) {
  if ( method.type == "MC" ) {
    # report start
    cat("\n")
    print(paste0(paste(rep("#", times = 20), collapse = ""), " ",
                 "RUNNING [RANDOM Ft - ASIAN OPTION] WITH ",
                 toupper(method2fullname(method.type)), " ",
                 paste(rep("#", times = 20), collapse = "")))
    cat("\n")
    print(paste("START TO RUN", toupper(method2fullname(method.type)), "AT",
                as.character(Sys.time()), sep = " "))
    cat("\n")
    # run inner layer simulation
    run_inner_layer(
      sim.params$contract.params, sim.params$inner.params, method.type,
      sim.params$Ft.mat, NULL,
      save.to.file = ifelse(method.type == "MC",
                            "../data/outer_asian_result_MC.rds",
                            paste0("../data/outer_asian_result_",
                                   method2filesuffix(method.type), "_",
                                   sim.params$outer.params$num.ref, ".rds")),
      ncpus = ncpus)
    saveRDS(sim.params,
            file = ifelse(method.type == "MC",
                          "../data/outer_asian_result_MC_sim_params.rds",
                          paste0("../data/outer_asian_result_",
                                 method2filesuffix(method.type), "_",
                                 sim.params$outer.params$num.ref,
                                 "_sim_params.rds")))
    # report end
    print(paste("END RUNNING", toupper(method2fullname(method.type)), "AT",
                as.character(Sys.time()), sep = " "))
    cat("\n")
  }
}
# run experiments with different reference points
for ( num.ref in c(5, 20, 50, 100) ) {
  sim.params$outer.params$num.ref <- num.ref
  sim.params[["Ft.ref.mat"]] <-
    t(apply(sim.params[["Ft.mat"]], MARGIN = 1, # apply by row (by each replication)
            FUN = function(row) as.numeric(
              quantile(row, probs =
                         seq_len(sim.params$outer.params$num.ref) /
                         ( sim.params$outer.params$num.ref + 1 ),
                       type = 1))))
  # report start
  cat("\n")
  print(paste0(paste(rep("#", times = 20), collapse = ""), " ",
               "RUNNING [RANDOM Ft - ASIAN OPTION] WITH ", num.ref, " REF PTS ",
               paste(rep("#", times = 20), collapse = "")))
  cat("\n")
  for ( method.type in sim.params$sim.methods ) {
    if ( method.type != "MC" ) {
      # report method name
      print(paste("START TO RUN", toupper(method2fullname(method.type)), "AT",
                  as.character(Sys.time()), sep = " "))
      cat("\n")
      # run inner layer simulation
      run_inner_layer(
        sim.params$contract.params, sim.params$inner.params, method.type,
        sim.params$Ft.mat, sim.params$Ft.ref.mat,
        save.to.file = ifelse(method.type == "MC",
                              "../data/outer_asian_result_MC.rds",
                              paste0("../data/outer_asian_result_",
                                     method2filesuffix(method.type), "_",
                                     sim.params$outer.params$num.ref, ".rds")),
        ncpus = ncpus)
      saveRDS(sim.params,
              file = ifelse(method.type == "MC",
                            "../data/outer_asian_result_MC_sim_params.rds",
                            paste0("../data/outer_asian_result_",
                                   method2filesuffix(method.type), "_",
                                   sim.params$outer.params$num.ref,
                                   "_sim_params.rds")))
      # report end
      print(paste("END RUNNING", toupper(method2fullname(method.type)), "AT",
                  as.character(Sys.time()), sep = " "))
      cat("\n")
    }
  }
}
################################################################################

################################################################################
#################### Sample Recycling Use Case Illustration ####################
################################################################################
##### Asian Option #####
# set up parameters
sim.params <- list(
  "contract.params" = contract.params.asian,
  "inner.params"    = inner.params.asian,
  "outer.params"    = outer.params.asian,
  "sim.methods"     = c("lambda_true", "lambda_naive")
)
sim.params$inner.params$N.sim <- 1
set.seed(428571) # set seed for reproducibility
Ft.RNG <- get_Ft_RNG(sim.params$outer.params, sim.params$contract.params)
sim.params[["Ft.ref.mat"]] <- matrix(Ft.RNG(n = sim.params$outer.params$num.ref),
                                     nrow = sim.params$inner.params$N.sim)
skip_RNG_streams(ncpus, # work with L'Ecuyer-CMRG RNG to ensure reproducibility
                 envir = globalenv())
sim.params[["Ft.mat"]] <- matrix(sort(c(
  Ft.RNG(n = sim.params$outer.params$N.out - sim.params$outer.params$num.ref),
  sim.params[["Ft.ref.mat"]])),
  nrow = sim.params$inner.params$N.sim)
skip_RNG_streams(ncpus, # work with L'Ecuyer-CMRG RNG to ensure reproducibility
                 envir = globalenv())
# run experiments with different methods
for ( method.type in sim.params$sim.methods ) {
  if ( method.type != "MC" ) {
    # report start
    cat("\n")
    print(paste0(paste(rep("#", times = 20), collapse = ""), " ",
                 "RUNNING [USE CASE - ASIAN OPTION] WITH ",
                 toupper(method2fullname(method.type)), " ",
                 paste(rep("#", times = 20), collapse = "")))
    cat("\n")
    print(paste("START TO RUN", toupper(method2fullname(method.type)), "AT",
                as.character(Sys.time()), sep = " "))
    cat("\n")
    # run inner layer simulation
    run_inner_layer(
      sim.params$contract.params, sim.params$inner.params, method.type,
      sim.params$Ft.mat, sim.params$Ft.ref.mat,
      save.to.file = ifelse(method.type == "MC",
                            "../data/use_case_asian_result_MC.rds",
                            paste0("../data/use_case_asian_result_",
                                   method2filesuffix(method.type), "_",
                                   sim.params$outer.params$num.ref, ".rds")),
      ncpus = ncpus)
    saveRDS(sim.params,
            file = ifelse(method.type == "MC",
                          "../data/use_case_asian_result_MC_sim_params.rds",
                          paste0("../data/use_case_asian_result_",
                                 method2filesuffix(method.type), "_",
                                 sim.params$outer.params$num.ref,
                                 "_sim_params.rds")))
    # report end
    print(paste("END RUNNING", toupper(method2fullname(method.type)), "AT",
                as.character(Sys.time()), sep = " "))
    cat("\n")
  }
}
################################################################################

# print session info
print(sessionInfo())
