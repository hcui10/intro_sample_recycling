## This file contains implementations for relevant algorithms accompanying 
## [Introducing Sample Recycling Method]
################################################################################

################################################################################
############## Simulation: Geometric Brownian Motion Sample Path ###############
################################################################################
# simulate one path
sim_path <- function(dt, tau, F0, mu, sigma) { 
  Wt <- cumsum(rnorm(as.integer(tau/dt), mean = 0, sd = sqrt(dt)))
  timeline <- seq(from = dt, to = tau, by = dt)
  return( c(F0, F0*exp( (mu-0.5*sigma^2)*timeline + sigma*Wt )) )
}
################################################################################

################################################################################
### Theoretical Price for European and Geometric Average-Rate Asian Options ####
################################################################################
# European Option
european_option_price_theory <- function(S0, K, sigma, r, d, tau, option.type) {
  d1 <- ( log( S0 / K) + ( r - d + 0.5 * sigma^2 ) * tau ) / sigma / sqrt(tau)
  d2 <- d1 - sigma * sqrt(tau)
  S.disc <- S0 * exp( - d * tau )
  K.disc <- K  * exp( - r * tau )
  if (option.type == "call") return( S.disc * pnorm(d1) - K.disc * pnorm(d2) )
  if (option.type == "put")  return( K.disc * pnorm(-d2) - S.disc * pnorm(-d1) )
}
# Geometric Continuous Average-Rate Asian Options
asian_option_price_theory <- function(S0, K, sigma, r, d, tau, option.type) {
  d_adj <- 0.5 * ( r + d + sigma^2 / 6 )
  sigma_adj <- sigma / sqrt(3)
  return( 
    european_option_price_theory(S0, K, sigma_adj, r, d_adj, tau, option.type) )
}
# Geometric Average-Rate Asian Options with Discrete Sample Steps
asian_option_price_theory_discrete <- function(S0, K, sigma, r, d, tau, 
                                               option.type, n, j, S.arr = NULL) {
  # incorporating step params: n discretely sampled steps
  if (n == Inf) {
    tau_mu <- tau / 2
    tau_sigma <- tau / 3
    B <- 1
    
  } else {
    h <- tau / n
    tau_mu <- ( 1 - j/n ) * ( tau - h * (n-j-1) / 2 ) 
    tau_sigma <- tau * ( 1 - j/n )^2 - 
      (n-j) * (n-j-1) * (4 * n - 4 * j + 1) / 6 / n / n * h
    B <- ifelse( j == 0, 1, prod(S.arr[1:j] / S0)^(1/n) )
  }
  # core pricing components
  A <- exp( - r * (tau - tau_mu) - d * tau_mu - 
              sigma^2 * (tau_mu - tau_sigma) * 0.5 ) * B
  d2 <- ( log( S0 / K ) + ( r - d - 0.5 * sigma^2 ) * tau_mu + log(B) ) / 
    sigma / sqrt(tau_sigma)
  d1 <- d2 + sigma * sqrt(tau_sigma)
  # option type: call / put 
  if (option.type == "call") {
    omega <- 1
  } else if (option.type == "put") { 
    omega <- -1 
  } else { 
    omega <- NULL
  }
  return( omega * S0 * A * pnorm( omega * d1 ) - 
            omega * K * exp( - r * tau ) * pnorm( omega * d2 ) )
}
################################################################################

################################################################################
######################## Vanilla Monte Carlo Simulation ########################
################################################################################
# Wrapper to Compute Average: Arithmetic and Geometric
avg <- function(arr, avg.method) {
  if (avg.method == "arithmetic") {
    return( mean(arr) )
  } else if (avg.method == "geometric") {
    return( exp( mean( log(arr[arr > 0]) ) ) )
  }
}
# Loss Path Functional 
eval_path_loss <- function(sample.path, loss.type, option.type, params) {
  # params: 
  # - European: K
  # - Asian: avg.target, avg.method, avg.idx, K (if average price)
  if (loss.type == "European") { 
    if (option.type == "call") { 
      return( max(sample.path[length(sample.path)] - params$K, 0) ) 
    } else if (option.type == "put") { 
      return( max(params$K - sample.path[length(sample.path)], 0) )
    }
  } else if (loss.type == "Asian") {
    if (params$avg.target == "price") {
      price <- avg(sample.path[params$avg.idx], avg.method = params$avg.method)
      strike <- params$K
    } else if (params$avg.target == "strike") {
      price <- sample.path[length(sample.path)]
      strike <- avg(sample.path[params$avg.idx], avg.method = params$avg.method)
    }
    if (option.type == "call") return( max(price - strike, 0) ) 
    if (option.type == "put")  return( max(strike - price, 0) ) 
  }
}
# Monte Carlo Option Pricing 
option_price_MC <- function(S0, K, sigma, r, d, tau, 
                            loss.type, option.type, 
                            num.MC, dt, avg.target, avg.method, avg.step, 
                            ncpus = 1, timeit = TRUE) {
  # simulate sample paths
  sim.time <- system.time({
    sample.paths <- 
      simplify2array( # each col is a path
        parallel::mclapply(seq(num.MC), function(x) sim_path(dt, tau, S0, r - d, sigma), 
                           mc.cores = ncpus, mc.allow.recursive = FALSE))
  })
  # evaluate losses 
  loss.eval.time <- system.time({
    # subset elements to take average on
    avg.idx <- seq(from = 1, to = nrow(sample.paths), by = avg.step)[-1]
    MC.losses <- as.numeric(simplify2array(parallel::mclapply(
      data.frame(sample.paths), 
      FUN = function(sample.path) 
        eval_path_loss(sample.path, loss.type = loss.type, option.type = option.type, 
                       params = list(avg.target = avg.target, 
                                     avg.method = avg.method, 
                                     avg.idx = avg.idx, 
                                     K = K)), 
      mc.cores = ncpus, mc.allow.recursive = FALSE))) * exp( - r * tau )
  })
  # return MC results
  if (timeit) {
    return( list( est = mean(MC.losses), 
                  disc.losses = MC.losses, 
                  samples = sample.paths, 
                  sim.time = sim.time, 
                  loss.eval.time = loss.eval.time) )
  } else {
    return( list( est = mean(MC.losses), 
                  disc.losses = MC.losses, 
                  samples = sample.paths ) )
  }
}
################################################################################

################################################################################
########################### Density Ratio Estimation ###########################
################################################################################
# Density Ratio for Geometric Brownian Motion
GBM_lambda <- function(F.tar, F.ref, sigma, r, d, dt) { # true density ratio
  # derived parameters
  F.tar.meanlog <- log(F.tar) + (r - d - 0.5 * sigma^2) * dt
  F.tar.sdlog   <- sigma * sqrt(dt)
  F.ref.meanlog <- log(F.ref) + (r - d - 0.5 * sigma^2) * dt
  F.ref.sdlog   <- F.tar.sdlog 
  # construct density ratio function 
  lambda_true <- function(x) 
    dlnorm(x, meanlog = F.tar.meanlog, sdlog = F.tar.sdlog) / 
    dlnorm(x, meanlog = F.ref.meanlog, sdlog = F.ref.sdlog)
  return( list("tar_meanlog" = F.tar.meanlog, "tar_sdlog" = F.tar.sdlog, 
               "ref_meanlog" = F.ref.meanlog, "ref_sdlog" = F.ref.sdlog, 
               "lambda" = lambda_true) )
}
# Density Ratio Estimation - Naive Stepwise
lambda_step_approx <- function(x_nu, x_de, n_block, x_min = -Inf, x_max = Inf) {
  # use denominator for breaks 
  breaks <-
    c(x_min, 
      as.numeric(head( # remove first and last (sample min and max) from quantiles 
        quantile(x_de, probs = seq(from = 0, to = 1,length.out = n_block + 1)[-1]), 
        -1)),
      x_max)
  # construct ratio
  n_nu <- hist(x_nu, breaks = breaks, plot = FALSE)$counts
  n_de <- hist(x_de, breaks = breaks, plot = FALSE)$counts
  ratio <- n_nu / n_de
  # Remove NAs, NaNs and Infs due to 0 counts
  ratio[is.na(ratio)] <- 0
  ratio[is.nan(ratio)] <- 0
  ratio[is.infinite(ratio)] <- 0
  # return estimated stepwise function 
  lambda <- approxfun(x = breaks, y = c(ratio, ratio[n_block]), 
                      # return NA for points outside the interval [min(x), max(x)]
                      rule = 1, 
                      # stepwise constant
                      method = "constant")
  return(lambda)
}
# Density Ratio Estimation Wrapper Function 
ratio_est <- function(classifier.type, x_nu, x_de, params = NULL) {
  # params: 
  # - GBM_true: F.tar, F.ref, sigma, r, d, dt
  # - naive_stepwise: n_block, x_min, x_max
  
  # fit data
  if ( classifier.type == "GBM_true" ) {
    # unpack params
    F.tar <- params$F.tar
    F.ref <- params$F.ref
    sigma <- params$sigma
    r <- params$r
    d <- params$d
    dt <- params$dt
    GBM_true <- GBM_lambda(F.tar, F.ref, sigma, r, d, dt)
    return( GBM_true$lambda )
  } else if ( classifier.type == "naive_stepwise" ) {
    naive_est <- 
      lambda_step_approx(x_nu, x_de, n_block = params$n_block, 
                         x_min = params$x_min, 
                         x_max = params$x_max)
    return( naive_est )
  } 
}
# map reference index to target indices (VERSION 1)
find_tar_idx <- function(Ft.ref.axis, Ft.axis) { # midpoint references
  # get distance to each ref pt
  dist.mat <- sapply(Ft.ref.axis, function(Ft.ref) abs(Ft.axis-Ft.ref))
  # get the ref point with the min distance
  dist.min.arr <- apply(dist.mat, MARGIN = 1, # apply by row
                        FUN = function(row) which.min(row)[1]) # pickfirst if tie
  # construct return list
  idx.tars.ls <- list()
  for (i in seq_len(length(Ft.ref.axis))) {
    idx.ref <- which(Ft.axis == Ft.ref.axis[i])
    idx.tars <- which(dist.min.arr == i)
    if ( idx.ref %in% idx.tars && length(idx.tars) == 1 ) {
      idx.tars.ls[[toString(idx.ref)]] <- NA
    } else {
      idx.tars.ls[[toString(idx.ref)]] <- idx.tars[!(idx.tars %in% c(idx.ref))]
    }
  }
  return(idx.tars.ls)
}
# # map reference index to target indices (VERSION 2)
# find_tar_idx <- function(Ft.ref.axis, Ft.axis) { # look left references
#   # construct return list
#   idx.tars.ls <- list()
#   for (i in seq_len(length(Ft.ref.axis))) {
#     idx.ref <- which(Ft.axis == Ft.ref.axis[i])
#     
#     if ( i == 1 ) tar.start <- 1
#     else tar.start <- which(Ft.axis == Ft.ref.axis[i-1]) + 1
#     
#     if ( tar.start + 1 == idx.ref ) idx.tars.ls[[toString(idx.ref)]] <- NA
#     else idx.tars.ls[[toString(idx.ref)]] <- as.integer(
#       seq(from = tar.start, to = idx.ref - 1, by = 1))
#   }
#   return(idx.tars.ls)
# }
# # map target index to reference index (VERSION 3)
# find_ref_idx <- function(Ft.ref.axis, Ft.axis) { # midpoint references
#   # get distance to each ref pt
#   dist.mat <- sapply(Ft.ref.axis, function(Ft.ref) abs(Ft.axis-Ft.ref))
#   # get the ref point with the min distance
#   dist.min.arr <- apply(dist.mat, MARGIN = 1, # apply by row
#                         FUN = function(row) which.min(row)[1]) # pickfirst if tie
#   names(dist.min.arr) <- as.character(seq_len(length(dist.min.arr)))
#   return(dist.min.arr)
# }
################################################################################

################################################################################
########################### Sample Recycling Method ############################
################################################################################
option_price_sample_recycle <- function(lambda.est, sample.paths.ref, disc.losses.ref) {
  F.test <- matrix(sample.paths.ref[2,], ncol = 1, byrow = TRUE)
  ratio.pred <- as.numeric(lambda.est(F.test))
  return( mean(ratio.pred * disc.losses.ref, na.rm = TRUE) )
}
################################################################################

################################################################################
######################### Helper Distribution Functions ########################
################################################################################
get_Ft_dist_params <- function(outer.params, contract.params) {
  Ft.dist.params <-           # theoretical Ft distribution params
    GBM_lambda(F.tar = outer.params$F0, 
               F.ref = outer.params$F0, 
               sigma = contract.params$sigma, 
               r  = outer.params$mu, 
               d  = contract.params$d, 
               dt = outer.params$t1)
  return(Ft.dist.params)
}
get_Ft_quantile <- function(outer.params, contract.params) {
  Ft.dist.params <- get_Ft_dist_params(outer.params, contract.params)
  Ft.quantile <- function(p)  # theoretical Ft quantiles
    qlnorm(p, meanlog = Ft.dist.params$tar_meanlog, sdlog = Ft.dist.params$tar_sdlog)
  return(Ft.quantile)
}
get_Ft_PDF <- function(outer.params, contract.params) {
  Ft.dist.params <- get_Ft_dist_params(outer.params, contract.params)
  Ft.PDF <- function(x)       # theoretical Ft PDF
    dlnorm(x, meanlog = Ft.dist.params$tar_meanlog, sdlog = Ft.dist.params$tar_sdlog)
  return(Ft.PDF)
}
get_Ft_CDF <- function(outer.params, contract.params) {
  Ft.dist.params <- get_Ft_dist_params(outer.params, contract.params)
  Ft.CDF <- function(q)       # theoretical Ft CDF
    plnorm(q, meanlog = Ft.dist.params$tar_meanlog, sdlog = Ft.dist.params$tar_sdlog)
  return(Ft.CDF)
}
get_Ft_RNG <- function(outer.params, contract.params) {
  Ft.dist.params <- get_Ft_dist_params(outer.params, contract.params)
  Ft.RNG <- function(n)       # random number generator from theoretical Ft distribution 
    rlnorm(n, meanlog = Ft.dist.params$tar_meanlog, sdlog = Ft.dist.params$tar_sdlog)
  return(Ft.RNG)
}
get_loss_func <- function(inner.params, contract.params) {
  return(switch(              # theoretical loss Lt
    contract.params$loss.type,
    "European" = function(Ft) european_option_price_theory(
      S0 = Ft, 
      K     = contract.params$K, 
      sigma = contract.params$sigma,
      r     = contract.params$r, 
      d     = contract.params$d,
      tau   = contract.params$tau, 
      option.type = contract.params$option.type),
    "Asian" = function(Ft) asian_option_price_theory_discrete(
      S0 = Ft, 
      K     = contract.params$K, 
      sigma = contract.params$sigma,
      r     = contract.params$r, 
      d     = contract.params$d,
      tau   = contract.params$tau, 
      option.type = contract.params$option.type,
      n = inner.params$n.dt, j = 0, S.arr = NULL)
  ))
}
get_Lt2Ft <- function(outer.params, inner.params, contract.params) {
  loss_func <- get_loss_func(inner.params, contract.params)
  Ft.quantile <- get_Ft_quantile(outer.params, contract.params)
  solve_Ft <- function(Lt.arr) sapply(Lt.arr, function(Lt) 
    uniroot(function(Ft) loss_func(Ft) - Lt, 
            interval = c(0, Ft.quantile(1-.Machine$double.eps^0.5)))$root )
}
get_Lt_PDF <- function(outer.params, inner.params, contract.params) {
  solve_Ft <- get_Lt2Ft(outer.params, inner.params, contract.params)
  Ft.PDF   <- get_Ft_PDF(outer.params, contract.params)
  Lt.PDF <- function(Lt.arr) { 
    if ( any(Lt.arr == 0) ) { # deal with inputs containing zeros 
      Lt.PDF.arr <- vector(mode = "numeric", length = length(Lt.arr))
      Lt.PDF.arr[Lt.arr == 0] <- 0
      if ( any(Lt.arr != 0) ) {
        Lt.PDF.arr[Lt.arr != 0] <- 
          Ft.PDF(solve_Ft(Lt.arr[Lt.arr != 0])) * 
          abs(numDeriv::grad(solve_Ft, x = Lt.arr[Lt.arr != 0]))
      }
      return(Lt.PDF.arr)
    } else {
      return( Ft.PDF(solve_Ft(Lt.arr)) * 
                abs(numDeriv::grad(solve_Ft, x = Lt.arr))  ) 
    }
  }
  return(Lt.PDF)
}
get_Lt_CDF <- function(outer.params, inner.params, contract.params) {
  solve_Ft <- get_Lt2Ft(outer.params, inner.params, contract.params)
  Ft.CDF   <- get_Ft_CDF(outer.params, contract.params)
  return( function(Lt.arr) Ft.CDF(solve_Ft(Lt.arr)) )
}
get_Lt_quantile <- function(outer.params, inner.params, contract.params) {
  loss_func <- get_loss_func(inner.params, contract.params)
  Ft.quantile <- get_Ft_quantile(outer.params, contract.params)
  return( function(probs) loss_func(Ft.quantile(probs)) )
}
################################################################################

################################################################################
################################# Risk Measures ################################
################################################################################
# Empirical Estimate of Risk Measures
risk_measure_est <- function(est.arr, risk.type, est.params = NULL) {
  # params: additional parameters for special risk types
  #   "Prob of Exceedance (POE)": thres K, array
  
  if ( risk.type == "VaR"  ) {
    n.axis <- min(101, length(est.arr))
    probs <- head(seq(from = 0, to = 1, length.out = n.axis)[-1], -1)
    return( as.numeric(quantile(est.arr, probs = probs)) )
  } else if ( risk.type == "CTE" ) { 
    n.axis <- min(101, length(est.arr))
    probs <- head(seq(from = 0, to = 1, length.out = n.axis)[-1], -1)
    VaRs <- as.numeric(quantile(est.arr, probs = probs))
    CTEs <- sapply(VaRs, function(VaR) 
      weighted.mean(est.arr, est.arr > VaR))
    return( CTEs )
  } else if ( risk.type == "CVaR" ) { 
    n.axis <- min(101, length(est.arr))
    probs <- head(seq(from = 0, to = 1, length.out = n.axis)[-1], -1)
    VaRs <- as.numeric(quantile(est.arr, probs = probs))
    CTEs <- sapply(VaRs, function(VaR) 
      weighted.mean(est.arr, est.arr > VaR))
    return( CTEs - VaRs )
  } else if ( risk.type == "POE" ) {
    return( as.numeric(sapply(est.params$K, function(K) 
      mean(est.arr > K))) )
  }
}
# Theoretical Functions of Risk Measures 
get_VaR_func <- function(outer.params, inner.params, contract.params) {
  Ft.quantile <- get_Ft_quantile(outer.params, contract.params)
  loss_func <- get_loss_func(inner.params, contract.params)
  return( function(probs) loss_func(Ft.quantile(probs)) )
}
get_CTE_func <- function(outer.params, inner.params, contract.params) {
  VaR_func <- get_VaR_func(outer.params, inner.params, contract.params)
  CTE_func <- function(probs) # theoretical conditional tail expectation (CTE)
    sapply(probs, function(prob) 
      ifelse(prob == 1, Inf, 
             integrate(VaR_func, lower = prob, upper = 1, 
                       rel.tol = .Machine$double.eps^0.5)$value / (1 - prob) ))
  return(CTE_func)
}
get_CVaR_func <- function(outer.params, inner.params, contract.params) {
  VaR_func <- get_VaR_func(outer.params, inner.params, contract.params)
  CTE_func <- get_CTE_func(outer.params, inner.params, contract.params)
  # theoretical conditional Value-at-Risk (CVaR)
  return( function(probs) CTE_func(probs) - VaR_func(probs) )
}
get_POE_func <- function(outer.params, inner.params, contract.params) {
  solve_Ft <- get_Lt2Ft(outer.params, inner.params, contract.params)
  Ft.CDF <- get_Ft_CDF(outer.params, contract.params)
  # theoretical probability of exceedance (POE)
  return( function(K.arr) return( 1 - Ft.CDF(solve_Ft(K.arr)) ) )
}
################################################################################
