#' Parallel variance calculation
#'
#' A helper function for `var_analytic` that calculations the variance elements
#' of the covariance matrix.
#'
#' @param fmm A "fastFMM" object
#' @param ... Additional arguments
#'
#' @return Calculated variance at location s
#' @export

var_parallel <- function(fmm, ...) {
  UseMethod("var_parallel")
}

#' Parallel variance calculation for non-concurrent models
#'
#' Calculate the variance of coefficients with index s.
#'
#' @param fmm A "fastFMM" object
#' @param mum Output from the massively univariate step
#' @param s Integer index for variance calculations
#' @param Z Matrix
#' @param RHat Smoothed coefficients
#' @param HHat Smoothed variance coefficients
#' @param id_list List of unique IDs
#' @param obs_ind List of vectors of rows corresponding to each subject
#' @param res_template Index template
#'
#' @return List of entries at index `s` for XTV inverse, the variance of
#' beta-tilde, and trimmed H-hat.
#'
#' @method var_parallel fastFMM
#' @export

# AX: separate issue: bug catching in parallel
# AX: report failed indices of calculation (e.g., due to non-invertible matrices)
# AX: E.g., message, consider removal of terminal points on functional domain

var_parallel.fastFMM <- function(
  fmm,
  mum,
  s,
  Z,
  RHat,
  HHat,
  id_list,
  obs_ind,
  res_template
) {
  # Setup parameters
  V_subj_inv <- c()
  randintercept <- mum$randintercept
  designmat <- mum$designmat
  p <- nrow(mum$betaTilde)
  L <- length(fmm$argvals)

  ## Invert each block matrix, then combine them
  if (!randintercept) {
    cov_trimmed <- eigenval_trim(
      matrix(c(HHat[, s], 0)[res_template], ncol(res_template))
    )
    # Returns a matrix
    HHat_trim_s <- cov_trimmed
  }

  # store XT * Vinv * X
  XTVinvX <- matrix(0, nrow = p, ncol = p)
  # store XT * Vinv * Z
  XTVinvZ_i <- vector(length = length(id_list), "list")

  ## iterate for each subject
  for (id in id_list) {
    subj_ind <- obs_ind[[as.character(id)]]
    subj_ind <- subj_ind

    # AX: Check for compatibility with concurrent case
    # Check dimensionality of Z[s][subj_ind, , drop = FALSE]
    if (randintercept) {
      Ji <- length(subj_ind)
      V_subj <- matrix(HHat[1, s], nrow = Ji, ncol = Ji) +
        diag(RHat[s], Ji)
    } else {
      V_subj <- Z[subj_ind, , drop = FALSE] %*%
        Matrix::tcrossprod(cov_trimmed, Z[subj_ind, , drop = FALSE]) +
        diag(RHat[s], length(subj_ind))
    }

    # Rfast provides a faster matrix inversion
    V_subj_inv <- as.matrix(Rfast::spdinv(V_subj))

    XTVinvX <- XTVinvX +
      crossprod(matrix(designmat[subj_ind,], ncol = p), V_subj_inv) %*%
      matrix(designmat[subj_ind,], ncol = p)


    if (randintercept) {
      XTVinvZ_i[[as.character(id)]] <- crossprod(
        matrix(designmat[subj_ind,], ncol = p),
        V_subj_inv
      ) %*% matrix(1, nrow = Ji, ncol = 1)
    } else {
      XTVinvZ_i[[as.character(id)]] <- crossprod(
        matrix(designmat[subj_ind, ], ncol = p),
        V_subj_inv
      ) %*% Z[subj_ind, , drop = FALSE]
    }
  }

  betaTilde_theo_var_s <- as.matrix(Rfast::spdinv(XTVinvX))

  return(
    list(
      XTVinv = XTVinvZ_i,
      betaTilde_theo_var = betaTilde_theo_var_s
    )
  )
}

#' Parallel variance calculation for concurrent models
#'
#' Calculate the variance of coefficients with index s.
#'
#' @param fmm A "fastFMM" object
#' @param mum Output from the massively univariate step
#' @param s Integer index for variance calculations
#' @param Z Matrix. This is only included to be consisntent with the
#' non-concurrent call and does nothing.
#' @param RHat Smoothed coefficients
#' @param HHat Smoothed variance coefficients
#' @param id_list List of unique IDs
#' @param obs_ind List of vectors of rows corresponding to each subject
#' @param res_template Index template
#'
#' @return List of entries at index `s` for XTV inverse, the variance of
#' beta-tilde, and trimmed H-hat.
#'
#' @method var_parallel fastFMMconc
#' @export

var_parallel.fastFMMconc <- function(
  fmm,
  mum,
  s,
  Z,
  RHat,
  HHat,
  id_list,
  obs_ind,
  res_template
) {
  # Setup parameters
  V_subj_inv <- c()
  randintercept <- mum$randintercept
  designmat <- mum$designmat[[s]]
  Z <- mum$ztlist[[s]]
  p <- nrow(mum$betaTilde)
  L <- length(fmm$argvals)

  ## Invert each block matrix, then combine them
  if (!randintercept) {
    cov_trimmed <- eigenval_trim(
      matrix(c(HHat[, s], 0)[res_template], ncol(res_template))
    )
    # Returns a matrix
    HHat_trim_s <- cov_trimmed
  }

  # store XT * Vinv * X
  XTVinvX <- matrix(0, nrow = p, ncol = p)
  # store XT * Vinv * Z
  XTVinvZ_i <- vector(length = length(id_list), "list")

  ## iterate for each subject
  for (id in id_list) {
    subj_ind <- obs_ind[[as.character(id)]]
    subj_ind <- subj_ind

    if (randintercept) {
      Ji <- length(subj_ind)
      V_subj <- matrix(HHat[1, s], nrow = Ji, ncol = Ji) +
        diag(RHat[s], Ji)
    } else {
      V_subj <- Z[subj_ind, , drop = FALSE] %*%
        Matrix::tcrossprod(cov_trimmed, Z[subj_ind, , drop = FALSE]) +
        diag(RHat[s], length(subj_ind))
    }

    # Rfast provides a faster matrix inversion
    V_subj_inv <- as.matrix(Rfast::spdinv(V_subj))

    XTVinvX <- XTVinvX +
      crossprod(matrix(designmat[subj_ind,], ncol = p), V_subj_inv) %*%
      matrix(designmat[subj_ind,], ncol = p)

    if (randintercept) {
      XTVinvZ_i[[as.character(id)]] <- crossprod(
        matrix(designmat[subj_ind,], ncol = p),
        V_subj_inv
      ) %*% matrix(1, nrow = Ji, ncol = 1)
    } else {
      XTVinvZ_i[[as.character(id)]] <- crossprod(
        matrix(designmat[subj_ind, ], ncol = p),
        V_subj_inv
      ) %*% Z[subj_ind, , drop = FALSE]
    }
  }

  betaTilde_theo_var_s <- as.matrix(Rfast::spdinv(XTVinvX))

  return(
    list(
      XTVinv = XTVinvZ_i,
      betaTilde_theo_var = betaTilde_theo_var_s
    )
  )
}

#' Parallelization wrapper
#'
#' To parallelize, the method needs to be dispatched properly. The normal
#' dispatch will pick the wrong argument, so we can use a wrapper. This
#' also has error handling capabilities to report the index of failure.
#'
#' @param s Integer index for variance calculations
#' @param fmm A "fastFMM" object
#' @param mum Output from the massively univariate step
#' @param Z Matrix. This is only included to be consisntent with the
#' non-concurrent call and does nothing.
#' @param RHat Smoothed coefficients
#' @param HHat Smoothed variance coefficients
#' @param id_list List of unique IDs
#' @param obs_ind List of vectors of rows corresponding to each subject
#' @param res_template Index template
#'
#' @return List of entries at index `s` for XTV inverse, the variance of
#' beta-tilde, and trimmed H-hat.
#'
#' @export

var_par_wrapper <- function(
    s, fmm, mum, Z, RHat, HHat, id_list, obs_ind, res_template
) {
  print("test")
  tryCatch(
    res <- var_parallel(fmm, mum, s, Z, RHat, HHat, id_list, obs_ind, res_template),
    warning = function(w) {
      warning("Warning in variance calculation at location ", s)
      print(w)
    },
    error = function(e) {
      stop("Variance calculation failed at location ", s)
      print(e)
    }
  )
}
