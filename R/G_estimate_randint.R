#' Special case of estimating covariance of random components G(s1, s2)
#'
#' Estimates the covariance matrix G for random intercepts that occurs at Step 3
#' of the FUI method. A helper function for `fui`.
#'
#' @param data A data frame containing all variables in formula
#' @param out_index Indices that contain the outcome variables
#' @param designmat Design matrix of the linear models
#' @param betaHat Estimated functional fixed effects
#' @param silent Whether to print the step description during calculations.
#' Defaults to `TRUE`.
#'
#' @return An estimation of the G matrix
#'
#' @importFrom Matrix crossprod
#'
#' @export

# TODO: check function dependencies

# Derive covariance estimates of random components: G(s1,s2)
### Create a function that estimates covariance G for random intercepts

G_estimate_randint <- function(
  data,
  out_index,
  designmat,
  betaHat,
  silent = TRUE
) {

  if(silent == FALSE)
    print("Step 3.1.1: Method of Moments Covariance Estimator Random Intercept")

  L <- length(out_index)
  GTilde <- matrix(NA, nrow = L, ncol = L)
  vdm <- crossprod(betaHat, var(designmat) %*% betaHat)
  d_temp <- data[,out_index]

  for(i in 1:L) {
    bhatVdm <- vdm[,i]
    d_temp_i <- d_temp[,i]
    res_temp <- GTilde[i,]
    for(j in 1:L){
      res_temp[j] <- stats::cov(d_temp_i, d_temp[,j], use = "pairwise.complete.obs") - bhatVdm[j]
    }
    GTilde[i,] <- res_temp
  }

  return(GTilde)
}
