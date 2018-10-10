## This file contains helper utility functions for other code accompanying 
## [Introducing Sample Recycling Method]
################################################################################

# from method to classifier type
method2classifier <- function(method.type) {
  return(
    switch(method.type, 
           "MC" = "MC", 
           "lambda_true" = "GBM_true", 
           "lambda_naive" = "naive_stepwise")
  )
}

# from method to file name suffix
method2filesuffix <- function(method.type) {
  return(
    switch(method.type, 
           "MC" = "MC", 
           "lambda_true" = "recy", 
           "lambda_naive" = "recy_naive")
  )
}

# from method to full method name
method2fullname <- function(method.type) {
  return(
    switch(method.type, 
           "MC" = "Monte Carlo Simulation", 
           "lambda_true" = "Sample Recycling Method", 
           "lambda_naive" = "Sample Recycling Method with Density Ratio Estimation")
  )
}

# from contract type to abbreviation
contract2abbr <- function(contract.type) {
  return(
    switch(contract.type, 
           "Asian" = "asian", 
           "European" = "euro")
  )
}

# summary statistics 
RMSE <- function(arr, true.val, na.rm = TRUE) {
  return( sqrt(mean( (arr-true.val)^2, na.rm = na.rm)) )
} 
MAPE <- function(arr, true.val, na.rm = TRUE) {
  return( median(abs(arr/true.val - 1), na.rm = na.rm) )
} 
Bias <- function(arr, true.val, na.rm = TRUE) {
  return( mean(arr, na.rm = na.rm) - true.val )
} 
summarize <- function(arr, true.val, na.rm = TRUE) {
  quartiles <- as.numeric(
    quantile(arr, probs = c(0, 0.25, 0.50, 0.75, 1)))
  summary_arr <- c(true.val, 
                   mean(arr, na.rm = na.rm),  
                   ifelse(is.na(true.val), NA, 
                          Bias(arr, true.val, na.rm = na.rm)), 
                   ifelse(is.na(true.val), NA, 
                          RMSE(arr, true.val, na.rm = na.rm) / true.val), 
                   ifelse(is.na(true.val), NA, 
                          MAPE(arr, true.val, na.rm = na.rm)), 
                   quartiles[1], quartiles[2], quartiles[3], 
                   quartiles[4], quartiles[5], 
                   as.integer(sum(is.na(arr))) )
  names(summary_arr) <- 
    c("True", "Mean", "Bias", "RMSE", "MAPE", 
      "Min", "1st Quartile", "Median", "3rd Quartile", "Max", 
      "Num of NAs")
  return(summary_arr)
}

# skip random number streams
skip_RNG_streams <- function(ncpus, envir) { # adpated from parallel package vignette
  s <- .Random.seed
  for (i in seq_len(ncpus)) s <- parallel::nextRNGStream(s)
  assign(".Random.seed", s, pos = envir)
}
