#' Estimate non-negative diagonal terms on G matrix
#'
#' Helper function for `G_estimate`. Uses least squares under non-negativity
#' constraints, mainly relying on `nnls` capability from `lsei`.
#'
#' @param data Data frame containing all predictor and outcome variables
#' @param out_index Indices of outcome variables in `data`
#' @param data_cov (unsure) Covariance of variables
#' @param RE_table Data frame containing random effects and interactions
#' @param idx_lst (unsure) Column indices of random effects
#' @param designmat (unsure)
#' @param betaHat Estimates of coefficients of random effects
#' @param GTilde Current `GTilde` estimate, created as an intermediate in `G_estimate`
#' @param non_neg (unsure), defaults to 0
#' @param silent Whether to print the step. Defaults to `TRUE`.

cov.nnls <- function(
  data,
  out_index,
  data_cov,
  RE_table,
  idx_lst,
  designmat,
  betaHat,
  GTilde,
  non_neg = 0,
  silent = TRUE
) {

  # Create helpful data frame of only outcomes
  d_temp <- data[, out_index]

  if (non_neg == 1) {

    if(silent == FALSE) print("Step 3.1.2: NNLS 1")

    # put constraints on EVERY coef corresponding to columns for one random effect
    ncol_Z <- ncol(data_cov$Z)
    var_term_idx <- which(is.na(RE_table$var2)) # find variance terms (that need non-negativity)

    # Initial indices
    idx_start_end <- rep(NA, ncol_Z)
    # concatenate all the terms that correspond to non-negative constraints
    non_negIdx <- do.call(c, lapply(var_term_idx, function(xx) idx_lst[[xx]]))
    # indices of unconstrained coefficients
    kk <- ncol_Z - length(non_negIdx)
    # put remaining terms after
    idx_start_end[1:kk] <- seq(1,ncol_Z)[-non_negIdx]
    # put unconstrained terms first
    idx_start_end[(kk+1):ncol_Z] <- non_negIdx
    # reorder to put covariates associated with non-negative coefs in first columns
    XX <- as.matrix(data_cov$Z[,idx_start_end])
    bHat <- rep(NA, dim(GTilde)[1])

    for (i in 1:L) {
      # Matrix multiply residuals
      YYi <- (d_temp[,i] - designmat %*% betaHat[,i])^2
      # NNLS step
      bHat[idx_start_end] <- lsei::pnnls(a = XX, b = YYi, k = kk)$x
      # average over coefficients corresponding to same random effect term
      GTilde[,i,i]  <- sapply(idx_lst,  function(x) mean(bHat[x]) )
    }

  } else if(non_neg == 2) {
    if(silent == FALSE) print("Step 3.1.2: NNLS 2")

    # put constraints on AVERAGE over coefs corresponding to columns for one random effect
    ncol_Z <- ncol(data_cov$Z)
    ff <- rep(0, nrow(RE_table)) # non-negativity vector
    eMat <- matrix(0, nrow = nrow(RE_table), ncol = ncol_Z) # constraint vector to be matrix below ( initially make all 0s so we do not place constraints on terms that can be negative)
    var_term_idx <- which(is.na(RE_table$var2)) # find variance terms (that need non-negativity)

    for (ii in var_term_idx) {
      eMat[ii, idx_lst[[ii]]] <- 1 # use these sum (and/or average) to enforce constraint on average
    }

    XX <- as.matrix(data_cov$Z) # reorder to put covariates associated with non-negative coefs in first columns and then transpose design mat for package

    for (i in 1:L) {
      YYi <- (d_temp[,i] - designmat %*% betaHat[,i])^2 # outcome of product of residuals
      bHat <- lsei::lsei(a = XX, b = YYi, e = eMat, f = ff) # nnls -- allows tiny negative values due to error
      GTilde[,i,i]  <- sapply(idx_lst,  function(x) mean(bHat[x]) ) # average over coefficients corresponding to same random effect term
    }

  }

  return(GTilde)
}
