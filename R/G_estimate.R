#' Generic "G_estimate" dispatch
#'
#' Estimates covariance of random components. The related function for the
#' random intercept flag (`G_estimate_randint`) does not require a generic
#' because the operations are the same for concurrent and non-concurrent models.
#'
#' @param fmm "An object that inherits from the "massmm" class"fastFMM object.
#' @param ... Additional arguments
#'
#' @return An estimation of the G matrix
#' @export

G_estimate <- function(fmm, ...) {
  UseMethod("G_estimate")
}

#' Estimate covariance of random components G(s1, s2), non-concurrent
#'
#' Estimates the covariance matrix G for random intercepts that occurs at Step 3
#' of the FUI method. Applies when `G_generate` cannot provide an analytic
#' solution.
#'
#' A helper function for `var_analytic`.
#'
#' @param fmm "fastFMM" object.
#' @param mum Result from massive univariate step
#' @param betaHat Estimated functional fixed effects
#' @param HHat Estimate H(s)
#' @param non_neg Integer control for non-negativity calculation
#' @param MoM Controls method of moments estimator
#' @param silent Whether to print the step description during calculations.
#' Defaults to `TRUE`.
#'
#' @method G_estimate fastFMM
#' @return An estimation of the G matrix
#'
#' @importFrom Matrix crossprod tcrossprod
#' @importFrom MASS ginv
#' @export

G_estimate.fastFMM <- function(
  fmm, mum, betaHat, HHat, non_neg, MoM, silent
) {
  if(silent == FALSE)
    print("Step 3.1.1: Method of Moments Covariance Estimator")

  data <- fmm$data
  L <- length(fmm$argvals)
  data_cov <- G_generate(fmm, mum, MoM)
  GTilde <- array(NA, dim = c(nrow(HHat), L, L))
  idx_lst <- data_cov$idx_lst
  out_index <- mum$out_index
  designmat <- mum$designmat
  ztlist <- mum$ztlist
  z_names <- names(ztlist)
  # this has concatenated design matrix and sums over columns for ID variables
  # do.call(cbind, ztlist)
  Z <- as.matrix(data_cov$Z_orig)

  # first part of OLS
  if (MoM == 2) {
    # First part of OLS expression
    XTXX <- as.matrix(
      Matrix::tcrossprod(
        MASS::ginv(
          as.matrix(Matrix::crossprod(data_cov$Z))
        ),
        data_cov$Z)
    )
  } else if (MoM == 1) {
    # function to join matrices
    mat_concat <- function(yy, xx) {
      if (length(xx) > 1)
        return(rowSums(yy[,xx]))
      else
        return(yy[,xx])
    }

    # sum across columns of Z associated with same random effect
    ZZ <- do.call(
      cbind,
      lapply(
        idx_lst,
        function(xx)
          mat_concat(yy=data_cov$Z, xx=xx)
      )
    )

    # first part of OLS expression
    # AX: MoM == 1 for concurrent case -> speedup?
    # AX: RFast may have options
    XTXX <- Matrix::tcrossprod( MASS::ginv( Matrix::crossprod(ZZ) ), ZZ)
    rm(ZZ)

  }

  # save design matrix for random effects
  idx_vec_zlst <- sapply(ztlist, ncol) # number of cols in each matrix
  idx_vec_zlst <- c(0, cumsum(idx_vec_zlst)) # vector of indices

  d_temp <- data[,out_index]
  if (MoM == 2) {
    for (i in 1:L) {
      YYi <- (d_temp[,i] - designmat %*% betaHat[, i])
      for (j in i:L) {
        # outcome of product of residuals
        YYj <- YYi * (d_temp[,j] - designmat %*% betaHat[, j])
        # coefficients from OLS with pseudo-inverse
        bHat <- XTXX %*% YYj
        GTilde[, i, j] <- GTilde[, j, i] <- sapply(
          idx_lst, function(x) mean(bHat[x])
        )
      }
    }
    rm(bHat)
  } else if (MoM == 1) {
    for (i in 1:L) {
      YYi <- (d_temp[,i] - designmat %*% betaHat[,i])
      for (j in i:L) {
        # outcome of product of residuals
        YYj <- YYi * (d_temp[,j] - designmat %*% betaHat[,j])
        # coefficients from OLS with pseudo-inverse
        GTilde[,i,j] <- GTilde[,j,i] <- XTXX %*% YYj
      }
    }
  }
  rm(YYi, YYj, XTXX)

  # non-negative least squares for estimation of non-diagonal terms
  if (non_neg != 0) {

    if (MoM == 1 & non_neg == 2) {
      message(paste(
        "Method of Moments approach 1 estimator can only use NNLS estimation",
        "scheme 1. Proceeding with NNLS-1"))
      non_neg <- 1
    }

    GTilde <- cov_nnls(
      mum = mum,
      data = data,
      L = L,
      data_cov = data_cov,
      betaHat = betaHat,
      GTilde = GTilde,
      non_neg = non_neg,
      silent = silent
    )
  }

  dimnames(GTilde)[[1]] <- rownames(HHat) # use names so lfosr_cov_organize() function below knows how to organize sub matrices based on names

  return(list(GTilde = GTilde, data_cov = data_cov))
}

#' Estimate covariance of random components G(s1, s2), non-concurrent
#'
#' Estimates the covariance matrix G for random intercepts that occurs at Step 3
#' of the FUI method. Applies when `G_generate` cannot provide an analytic
#' solution.
#'
#' A helper function for `var_analytic`.
#'
#' @param fmm "fastFMM" object.
#' @param mum Result from massive univariate step
#' @param betaHat Estimated functional fixed effects
#' @param HHat Estimate H(s)
#' @param non_neg Integer control for non-negativity calculation
#' @param MoM Controls method of moments estimator
#' @param silent Whether to print the step description during calculations.
#' Defaults to `TRUE`.
#'
#' @method G_estimate fastFMMconc
#' @return An estimation of the G matrix
#'
#' @importFrom Matrix crossprod tcrossprod
#' @importFrom MASS ginv
#' @export

G_estimate.fastFMMconc <- function(
  fmm, mum, betaHat, HHat, non_neg, MoM, silent
) {
  if(silent == FALSE)
    print("Step 3.1.1: Method of Moments Covariance Estimator")

  # Dummy
  data <- fmm$data
  data_cov <- G_generate(fmm, mum, 1, 1)
  L <- length(fmm$argvals)
  GTilde <- array(NA, dim = c(nrow(HHat), L, L))
  idx_lst <- data_cov$idx_lst
  out_index <- fmm$out_index
  designmat <- mum$designmat
  ztlist <- mum$ztlist
  z_names <- names(ztlist)
  # this has concatenated design matrix and sums over columns for ID variables
  # do.call(cbind, ztlist)
  # Z <- as.matrix(data_cov$Z_orig)

  d_temp <- data[, out_index]

  # save design matrix for random effects
  idx_vec_zlst <- sapply(ztlist, ncol) # number of cols in each matrix
  idx_vec_zlst <- c(0, cumsum(idx_vec_zlst)) # vector of indices

  for (i in 1:L) {
    YYi <- (d_temp[, i] - designmat[[i]] %*% betaHat[, i])

    for(j in i:L){
      YYj <- YYi * (d_temp[, j] - designmat[[j]] %*% betaHat[, j]) # outcome of product of residuals
      data_cov <- G_generate(fmm, mum, i, j)
      XTXX <- as.matrix(tcrossprod(MASS::ginv(as.matrix(crossprod(data_cov$Z))), data_cov$Z))

      bHat <- XTXX %*% YYj  # coefficients from OLS with pseudo-inverse
      GTilde[, i, j] <- GTilde[, j, i] <- sapply(
        idx_lst, function(x) mean(bHat[x])
      )
    }
  }

  rm(bHat)
  rm(YYi, YYj, XTXX)

  # non-negative least squares for estimation of non-diagonal terms
  # AX: I'm actually not sure if this will work with the concurrent output
  # if (non_neg != 0) {
  #
  #   if (MoM == 1 & non_neg == 2) {
  #     message(paste(
  #       "Method of Moments approach 1 estimator can only use NNLS estimation",
  #       "scheme 1. Proceeding with NNLS-1"))
  #     non_neg <- 1
  #   }
  #
  #   GTilde <- cov_nnls(
  #     mum = mum,
  #     data = data,
  #     L = L,
  #     data_cov = data_cov,
  #     betaHat = betaHat,
  #     GTilde = GTilde,
  #     non_neg = non_neg,
  #     silent = silent
  #   )
  # }

  # use names so lfosr_cov_organize() function below knows how to organize
  # sub matrices based on names
  dimnames(GTilde)[[1]] <- rownames(HHat)

  return(list(GTilde = GTilde, data_cov = data_cov))
}


#' Estimate non-negative diagonal terms on G matrix
#'
#' Helper function for `G_estimate`. Uses least squares under non-negativity
#' constraints, mainly relying on `nnls` capability from `lsei`.
#'
#' @param mum Output of massively univariate step
#' @param data Data frame containing all predictor and outcome variables
#' @param L Integer dimension of functional domain
#' @param data_cov Output of `G_generate`
#' @param betaHat Estimates of coefficients of random effects
#' @param GTilde Current `GTilde` estimate, created as an intermediate.
#' @param non_neg Integer control for non-negativity step, defaults to 0
#' @param silent Whether to print the step. Defaults to `TRUE`.
#'
#' @return A new GTilde matrix with enforced non-negativity
#'
#' @importFrom lsei pnnls

cov_nnls <- function(
  mum,
  data,
  L,
  data_cov,
  betaHat,
  GTilde,
  non_neg = 0,
  silent = TRUE
) {

  out_index <- mum$out_index
  RE_table <- mum$varcorr_df
  idx_lst <- data_cov$idx_lst
  designmat <- mum$designmat
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
    XX <- as.matrix(data_cov$Z[, idx_start_end])
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
      YYi <- (d_temp[, i] - designmat %*% betaHat[, i])^2 # outcome of product of residuals
      bHat <- lsei::lsei(a = XX, b = YYi, e = eMat, f = ff) # nnls -- allows tiny negative values due to error
      GTilde[, i, i]  <- sapply(idx_lst,  function(x) mean(bHat[x]) ) # average over coefficients corresponding to same random effect term
    }

  }

  return(list(GTilde = GTilde, data_cov = data_cov))
}
