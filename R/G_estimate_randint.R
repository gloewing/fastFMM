#' Special case of estimating covariance of random components G(s1, s2)
#'
#' Estimates the covariance matrix G at Step 3. If the only random effect is an
#' intercept, we can use a speedup.
#'
#' @param fmm "fastFMM" object.
#' @param designmat Design matrix of the linear models. A list if concurrent.
#' @param betaHat Estimated functional fixed effects
#' @param silent Whether to print the step description during calculations.
#' Defaults to `TRUE`.
#'
#' @return An estimation of the G matrix
#'
#' @importFrom Matrix crossprod
#' @importFrom stats cov var
#' @export

# Derive covariance estimates of random components: G(s1,s2)
### Create a function that estimates covariance G for random intercepts

G_estimate_randint <- function(
  fmm,
  designmat,
  betaHat,
  silent = TRUE
) {

  if(silent == FALSE)
    print("Step 3.1: MoM Covariance Estimator: Only Random Intercept Case")

  data <- fmm$data
  out_index <- fmm$out_index
  L <- length(fmm$argvals)
  # Check if concurrent
  if (!is.null(fmm$fun_covariates))
      # Should be identical if only random intercept
      designmat <- designmat[[1]]
  GTilde <- matrix(NA, nrow = L, ncol = L)
  vdm <- Matrix::crossprod(betaHat, stats::var(designmat) %*% betaHat)
  d_temp <- data[, out_index]

  for(i in 1:L) {
    bhatVdm <- vdm[,i]
    d_temp_i <- d_temp[,i]
    res_temp <- GTilde[i,]
    for(j in 1:L){
      res_temp[j] <- stats::cov(
        d_temp_i,
        d_temp[,j],
        use = "pairwise.complete.obs"
      ) - bhatVdm[j]
    }
    GTilde[i, ] <- res_temp
  }

  return(GTilde)
}
