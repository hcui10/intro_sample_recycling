## This file contains helper functions to generate tables accompanying 
## [Introducing Sample Recycling Method]
################################################################################

# data validation 
check_data <- function(sim.results, F.axis, F.ref.axis) {
  # get partition indices
  idx.ref <- which(!is.na(match(F.axis, F.ref.axis)))
  idx.tar <- which(is.na(match(F.axis, F.ref.axis)))
  sim.methods <- names(sim.results)
  # loop through all methods in the simulation result
  for (method.type in sim.methods) {
    # error flag and counter
    err.flag <- FALSE
    err.counter <- 0
    # VALUE: loss estimate
    if ( sum(!is.finite(sim.results[[method.type]][["loss_est"]])) != 0 ) {
      print(paste(method.type, "has invalid (Inf/NA) loss est values", sep = " "))
      err.flag <- TRUE
      err.counter <- err.counter + 1
    }
    # VALUE: sample average of estimated lambda
    if ( method.type != "MC") {
      if ( sum(!is.na(sim.results[[method.type]][["exp_lambda"]][,idx.ref])) != 0 ) {
        print(paste(method.type, "has non-NA values for exp_lambda on ref pts", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
      if ( sum(!is.finite(sim.results[[method.type]][["exp_lambda"]][,idx.tar])) != 0 ) {
        print(paste(method.type, "has invalid (Inf/NA) values for exp_lambda on ref pts", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
    }
    # TIME: inner loop simulation 
    if ( method.type == "MC" ) {
      if ( sum(!is.finite(sim.results[[method.type]][["time"]][["sim"]])) != 0 ) {
        print(paste(method.type, "has invalid (Inf/NA) sim time", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
    } else {
      if ( sum(!is.finite(sim.results[[method.type]][["time"]][["sim"]][,idx.ref])) != 0 ) {
        print(paste(method.type, "has invalid (Inf/NA) sim time on ref pts", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
      if ( sum(sim.results[[method.type]][["time"]][["sim"]][,idx.tar], na.rm = FALSE) != 0 ) {
        print(paste(method.type, "has invalid (Inf/NA/non-zero) sim time on tar pts", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
    }
    # TIME: loss evaluation
    if ( sum(!is.finite(sim.results[[method.type]][["time"]][["eval"]])) != 0 ) {
      print(paste(method.type, "has invalid (Inf/NA) loss eval time", sep = " "))
      err.flag <- TRUE
      err.counter <- err.counter + 1
    }
    # TIME: density ratio - simulation 
    if (method.type != "MC") {
      if (method.type == "lambda_true") {
        if ( sum(sim.results[[method.type]][["time"]][["ratio_sim"]], na.rm = FALSE) != 0 ) {
          print(paste(method.type, "has invalid (Inf/NA/non-zero) ratio_sim time", sep = " "))
          err.flag <- TRUE
          err.counter <- err.counter + 1
        }
      } else {
        if ( sum(!is.finite(sim.results[[method.type]][["time"]][["ratio_sim"]])) != 0 ) {
          print(paste(method.type, "has invalid (Inf/NA) ratio_sim time", sep = " "))
          err.flag <- TRUE
          err.counter <- err.counter + 1
        }
      }
    }
    # TIME: density ratio - estimation 
    if (method.type != "MC") {
      if ( sum(sim.results[[method.type]][["time"]][["ratio_est"]][,idx.ref], na.rm = FALSE) != 0 ) {
        print(paste(method.type, "has invalid (Inf/NA/non-zero) ratio_est time on ref pts", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
      if ( sum(!is.finite(sim.results[[method.type]][["time"]][["ratio_est"]][,idx.tar])) != 0 ) {
        print(paste(method.type, "has invalid (Inf/NA) ratio_est time on tar pts", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
    }
    # TIME: density ratio - evaluation 
    if (method.type != "MC") {
      if ( sum(sim.results[[method.type]][["time"]][["ratio_eval"]][,idx.ref], na.rm = FALSE) != 0 ) {
        print(paste(method.type, "has invalid (Inf/NA/non-zero) ratio_eval time on ref pts", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
      if ( sum(!is.finite(sim.results[[method.type]][["time"]][["ratio_val"]][,idx.tar])) != 0 ) {
        print(paste(method.type, "has invalid (Inf/NA) ratio_eval time on tar pts", sep = " "))
        err.flag <- TRUE
        err.counter <- err.counter + 1
      }
    }
    # output data validation summary 
    if (!err.flag) {
      print(paste(method.type, "--- ALL DATA VALIDATION CLEARED ! ---", sep = " "))
    } else {
      print(paste("--- DATA VALIDATION END :", err.counter, "ERROR(S) FOUND ---", sep = " "))
    }
  }
}

# generate time profile table (as a matrix)
make_time_table <- function(method.type, sim.results, F.axis, F.ref.axis) {
  if (method.type == "MC") {
    # get Monte Carlo time profile
    time.profile.mat <- t(sapply(
      list(sim.results[[method.type]][["time"]][["sim"]],    # inner loop sample path simulation 
           sim.results[[method.type]][["time"]][["eval"]],   # loss evaluation on sample paths
           sim.results[[method.type]][["time"]][["eval"]] +  # total time per L(F) estimation 
             sim.results[[method.type]][["time"]][["sim"]]), 
      FUN = function(time.mat) 
        c(mean(as.numeric(time.mat), na.rm = TRUE), 
          sqrt(var(as.numeric(time.mat), na.rm = TRUE)))
    ))
    # fill in names
    colnames(time.profile.mat) <- c("Mean", "SD")
    rownames(time.profile.mat) <- 
      c("Inner Loop Simulation", "Loss Evaluation", "Overall Average")
  } else {
    # partition indices
    idx.ref <- which(!is.na(match(F.axis, F.ref.axis)))
    idx.tar <- which(is.na(match(F.axis, F.ref.axis)))
    # get sample recycling time profile
    time.profile.mat <- t(sapply(
      list(
        # ref - inner loop sample path simulation 
        sim.results[[method.type]][["time"]][["sim"]][,idx.ref], 
        # ref - loss evaluation on sample paths 
        sim.results[[method.type]][["time"]][["eval"]][,idx.ref], 
        # ref - denominator sample generation for density ratio estimation
        sim.results[[method.type]][["time"]][["ratio_sim"]][,idx.ref], 
        # ref - total time total time per L(F) estimation 
        sim.results[[method.type]][["time"]][["sim"]][,idx.ref] + 
          sim.results[[method.type]][["time"]][["eval"]][,idx.ref] + 
          sim.results[[method.type]][["time"]][["ratio_sim"]][,idx.ref] + 
          sim.results[[method.type]][["time"]][["ratio_est"]][,idx.ref] + 
          sim.results[[method.type]][["time"]][["ratio_eval"]][,idx.ref], 
        # tar - numerator sample generation for density ratio estimation
        sim.results[[method.type]][["time"]][["ratio_sim"]][,idx.tar], 
        # tar - estimate density ratio function given tar (and ref) 
        sim.results[[method.type]][["time"]][["ratio_est"]][,idx.tar], 
        # tar - apply estimated density ratio function (given tar) on ref samples
        sim.results[[method.type]][["time"]][["ratio_eval"]][,idx.tar], 
        # tar - loss evaluation on tar by recycling ref losses 
        sim.results[[method.type]][["time"]][["eval"]][,idx.tar], 
        # tar - total time total time per L(F) estimation 
        sim.results[[method.type]][["time"]][["sim"]][,idx.tar] + 
          sim.results[[method.type]][["time"]][["eval"]][,idx.tar] + 
          sim.results[[method.type]][["time"]][["ratio_sim"]][,idx.tar] + 
          sim.results[[method.type]][["time"]][["ratio_est"]][,idx.tar] + 
          sim.results[[method.type]][["time"]][["ratio_eval"]][,idx.tar], 
        # OVERALL - average time oer L(F) estimation replicated N.sim times 
        rowMeans(sim.results[[method.type]][["time"]][["sim"]] + 
                   sim.results[[method.type]][["time"]][["eval"]] + 
                   sim.results[[method.type]][["time"]][["ratio_sim"]] + 
                   sim.results[[method.type]][["time"]][["ratio_est"]] + 
                   sim.results[[method.type]][["time"]][["ratio_eval"]])
      ), FUN = function(time.mat) 
        c(mean(as.numeric(time.mat), na.rm = TRUE), 
          sqrt(var(as.numeric(time.mat), na.rm = TRUE)))
    ))
    # special handle for true density ratio
    if (method.type == "lambda_true") {
      # no sample simulation for density ratio estimation 
      time.profile.mat[3,] <- time.profile.mat[5,] <- NA
    }
    # fill in names
    colnames(time.profile.mat) <- c("Mean", "SD")
    rownames(time.profile.mat) <- 
      c("Ref - Inner Loop Simulation", 
        "Ref - Loss Evaluation", 
        "Ref - Sample Simulation for Density Ratio Est", 
        "Ref - Partial Average", 
        "Tar - Sample Simulation for Density Ratio Est", 
        "Tar - Density Ratio Estimation", 
        "Tar - Density Ratio Evaluation", 
        "Tar - Loss Evaluation", 
        "Tar - Partial Average", 
        "Overall Average")
  }
  return(time.profile.mat)
}
