#' Fast Univariate Inference for Longitudinal Functional Models
#'
#' Fit a function-on-scalar regression model for longitudinal
#' functional outcomes and scalar predictors using the Fast Univariate
#' Inference (FUI) approach (Cui et al. 2022).
#'
#' The FUI approach comprises of three steps:
#' \enumerate{
#' \item Fit a univariate mixed model at each location of the functional domain,
#' and obtain raw estimates from massive models;
#' \item Smooth the raw estimates along the functional domain;
#' \item Obtain the pointwise and joint confidence bands using an analytic
#' approach for Gaussian data or Bootstrap for general distributions.
#' }
#'
#' For more information on each step, please refer to the FUI paper
#' by Cui et al. (2022).
#'
#' @param formula Two-sided formula object in lme4 formula syntax.
#' The difference is that the response need to be specified as a matrix
#' instead of a vector. Each column of the matrix represents one location
#' of the longitudinal functional observations on the domain.
#' @param data A data frame containing all variables in formula
#' @param family GLM family of the response. Defaults to \code{gaussian}.
#' @param var Logical, indicating whether to calculate and return variance
#' of the coefficient estimates. Defaults to \code{TRUE}.
#' @param analytic Logical, indicating whether to use the analytic inferenc
#' approach or bootstrap. Defaults to \code{TRUE}.
#' @param parallel Logical, indicating whether to do parallel computing.
#' Defaults to \code{FALSE}.
#' @param silent Logical, indicating whether to show descriptions of each step.
#' Defaults to \code{FALSE}.
#' @param argvals A vector containing locations of observations on the
#' functional domain. If not specified, a regular grid across the range of
#' the domain is assumed. Currently only supported for bootstrap
#' (\code{analytic=FALSE}).
#' @param nknots_min Minimal number of knots in the penalized smoothing for the
#' regression coefficients.
#' Defaults to \code{NULL}, which then uses L/2 where L is the dimension of the
#' functional domain.
#' @param nknots_min_cov Minimal number of knots in the penalized smoothing for
#' the covariance matrices.
#' Defaults to \code{35}.
#' @param smooth_method How to select smoothing parameter in step 2. Defaults to
#'  \code{"GCV.Cp"}
#' @param splines Spline type used for penalized splines smoothing. We use the
#' same syntax as the mgcv package. Defaults to \code{"tp"}.
#' @param design_mat Logical, indicating whether to return the design matrix.
#' Defaults to \code{FALSE}
#' @param residuals Logical, indicating whether to save residuals from
#' unsmoothed LME. Defaults to \code{FALSE}.
#' @param num_boots Number of samples when using bootstrap inference. Defaults
#' to 500.
#' @param boot_type Bootstrap type (character): "cluster", "case", "wild",
#' "reb", "residual", "parametric", "semiparametric". \code{NULL} defaults to
#' "cluster" for non-gaussian responses and "wild" for gaussian responses. For
#' small cluster (n<=10) gaussian responses, defaults to "reb".
#' @param seed Numeric value used to make sure bootstrap replicate (draws) are
#' correlated across functional domains for certain bootstrap approach
#' @param subj_id Name of the variable that contains subject ID.
#' @param num_cores Number of cores for parallelization. Defaults to 1.
#' @param caic Logical, indicating whether to calculate cAIC. Defaults to
#' \code{FALSE}.
#' @param REs Logical, indicating whether to return random effect estimates.
#' Defaults to \code{FALSE}.
#' @param non_neg 0 - no non-negativity constrains, 1 - non-negativity
#' constraints on every coefficient for variance, 2 - non-negativity on
#' average of coefficents for 1 variance term. Defaults to 0.
#' @param MoM Method of moments estimator. Defaults to 1.
#' @param impute_outcome Logical, indicating whether to impute missing outcome
#' values with FPCA. Defaults to \code{FALSE}.
#' @param concurrent Logical, indicating whether to fit a concurrent model.
#' Defaults to \code{FALSE}.
#'
#' @return A list containing:
#' \item{betaHat}{Estimated functional fixed effects}
#' \item{argvals}{Location of the observations}
#' \item{betaHat.var}{Variance estimates of the functional fixed effects
#' (if specified)}
#' \item{qn}{critical values used to construct joint CI}
#' \item{...}{...}
#'
#' @author Erjia Cui \email{ecui@@umn.edu}, Gabriel Loewinger
#' \email{gloewinger@@gmail.com}
#'
#' @references Cui, E., Leroux, A., Smirnova, E., Crainiceanu, C. (2022). Fast
#' Univariate Inference for Longitudinal Functional Models. \emph{Journal of
#' Computational and Graphical Statistics}, 31(1), 219-230.
#'
#' @export
#'
#' @import lme4
#' @import parallel
#' @import magrittr
#' @import dplyr
#' @import stringr
#' @import mgcv
#' @import refund
#' @importFrom MASS ginv
#' @import cAIC4
#' @importFrom lsei pnnls lsei
#' @import Matrix
#' @importFrom mvtnorm rmvnorm
#' @importFrom Rfast rowMaxs spdinv
#' @import progress
#'
#' @examples
#' library(refund)
#'
#' ## random intercept only
#' set.seed(1)
#' DTI_use <- DTI[DTI$ID %in% sample(DTI$ID, 10),]
#' fit_dti <- fui(
#'   cca ~ case + visit + sex + (1 | ID),
#'   data = DTI_use
#' )

fui <- function(
  formula,
  data,
  family = "gaussian",
  var = TRUE,
  analytic = TRUE,
  parallel = FALSE,
  silent = FALSE,
  argvals = NULL,
  nknots_min = NULL,
  nknots_min_cov = 35,
  smooth_method = "GCV.Cp",
  splines = "tp",
  design_mat = FALSE,
  residuals = FALSE,
  num_boots = 500,
  boot_type = NULL,
  seed = 1,
  subj_id = NULL,
  num_cores = 1,
  caic = FALSE,
  REs = FALSE,
  non_neg = 0,
  MoM = 1,
  concurrent = FALSE,
) {

  # 0. Setup ###################################################################

  # 0.1 Set parallelization ==================================================

  # If doing parallel computing, set up the number of cores
  if (parallel & !is.integer(num_cores)) {
    num_cores <- as.integer(round(parallel::detectCores() * 0.75))
    if (num_cores < 2)
      warning("Only 1 core detected for parallelization.")
    if (!silent)
      message(paste("No. cores used for parallelization:", num_cores))
  }

  # For non-Gaussian family, only do bootstrap inference
  if (family != "gaussian")
    analytic <- FALSE

  # 0.2  Organize the input from the model formula =============================

  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)

  # Stop if there are column names with "." to avoid issues with covariance
  # G() and H() calculations
  dep_str <- deparse(model_formula[3])
  if (grepl(".", dep_str, fixed = TRUE)) {
    # make sure it isn't just a call to all covariates with "Y ~. "
    # remove first character of parsed formula string and check
    dep_str_rm <- substr(dep_str, 3, nchar(dep_str))
    if (grepl(".", dep_str_rm, fixed = TRUE)) {
      stop(
        paste0(
          'Remove the character "." from all non-functional covariate names ',
          'and rerun fui()', '\n',
          '- E.g., change "X.1" to "X_1"', '\n',
          '- The string "." *should* be kept in functional outcome names ',
          '(e.g., "Y.1" *is* proper naming).'
        )
      )
    }
  }
  rm(dep_str)

  # Obtain the dimension of the functional domain
  # Find indices that start with the outcome name
  out_index <- grep(paste0("^", model_formula[2]), names(data))

  # Read the full list of columns or the matrix column
  # Set L depending on whether out_index is multiple columns or a matrix
  if (length(out_index) != 1) {
    # Multiple columns
    L <- length(out_index)
  } else {
    # Matrix column
    L <- ncol(data[, out_index])
  }

  # 0.3 Impute missing values ==================================================

  # Fill in missing values of functional outcome using FPCA
  # rows with missing outcome values

  # AX: Consider shortening this through writing a new function.

  missing_rows <- which(rowSums(is.na(data[, out_index])) != 0)
  if (length(missing_rows) != 0) {
    if(analytic & impute_outcome){
      message(
        paste(
          "Imputing", length(out_index),
          "values in functional response with longitudinal functional PCA"
        )
      )
      if (length(out_index) != 1) {
        tmp <- as.matrix(data[, out_index])
        tmp[which(is.na(tmp))] <- suppressWarnings(
          refund::fpca.face(
            tmp,
            argvals = argvals,
            knots = nknots_fpca
          )$Yhat[which(is.na(tmp))]
        )
        data[,out_index] <- tmp
      } else {
        data[, out_index][which(is.na(data[, out_index]))] <- suppressWarnings(
          refund::fpca.face(
            as.matrix(data[, out_index]),
            argvals = argvals,
            knots = nknots_fpca
          )$Yhat[which(is.na(data[, out_index]))]
        )
      }
    } else if (analytic & !impute_outcome) {
      message(
        paste(
          "Removing", missing_rows,
          "rows with missing functional outcome values.", "\n",
          "To impute missing outcome values with FPCA, set fui() argument:",
          "impute_outcome = TRUE"
        )
      )

      # remove data with missing rows
      data[, out_index] <- data[-missing_rows, out_index]
    }
  }

  # 0.4 Functional covariates ==================================================

  # Detect functional covariates
  x_names <- all.vars(formula)[-1]
  x_classes <- sapply(x_names, function(x) class(data[[x]]))
  x_ncols <- sapply(
    x_names,
    function(x) length(grep(paste0("^", x), names(data)))
  )
  func_covs <- NULL # Redundacy for clarity
  func_covs <- unique(c(x_names[x_classes == "AsIs"], x_names[x_ncols == L]))

  # Check for inconsistencies between concurrence arg and functional covariates
  if (concurrent & length(func_covs < 1)) {
    stop(
      "No functional covariates found for concurrent model fitting.", "\n",
      "Incoporate functional covariates or set `concurrent = F` in args."
    )
  } else if (!concurrent & length(func_covs) > 0) {
    warning(
      "Functional covariates detected: ",
      paste0(func_covs, collapse = ", "), "\n",
      "Execution will continue, but consider setting `concurrent = TRUE`."
    )
  }

  # 1. Massively univariate mixed models #######################################

  if (!silent)
    print("Step 1: Fit Massively Univariate Mixed Models")

  # 1.0 Model sanity checks ====================================================

  if (analytic & !is.null(argvals) & var)
    message(
      paste(
        "'argvals' argument is currently only supported for bootstrap.",
        "`argvals' ignored: fitting model to ALL points on functional domain"
      )
    )

  if (is.null(argvals) | analytic) {
    argvals <- 1:L
  } else {
    if (max(argvals) > L)
      stop(
        paste(
          "Maximum index specified in argvals is greater than",
          "the total number of columns for the functional outcome"
        )
      )
    L <- length(argvals)
  }

  if (family == "gaussian" & analytic & L > 400 & var)
    message(
      paste(
        "Your functional data is dense!",
        "Consider subsampling along the functional domain",
        "(i.e., reduce columns of outcome matrix)",
        "or using bootstrap inference."
      )
    )

  # Create a matrix to store AICs
  AIC_mat <- matrix(NA, nrow = L, ncol = 2)

  # 1.1 Fitting ================================================================

  # Initialize parameters within the univariate model
  # Create an object with the appropriate class: "unimm" or "unimm_conc"
  uni_model <- new_unimm(
    formula = formula,
    family = family,
    residuals = residuals,
    caic = caic,
    REs = REs,
    analytic = analytic,
    concurrent = concurrent,
    func_covs = func_covs
  )

  # Create a list of univariate models ("massively univariate")
  mum <- massmm(uni_model, argvals, data, parallel)

  # 2. Smoothing ###############################################################

  if (!silent) print("Step 2: Smoothing")

  # Penalized splines smoothing and extract components (analytic)
  # Number of knots for regression coefficients
  nknots <- min(round(L / 2), nknots_min)
  # Number of knots for covariance matrix
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- min(round(L / 2), 35)

  # Smoothing parameter, spline basis, penalty matrix (analytic)
  if (analytic) {
    p <- nrow(mum$betaTilde) # Number of fixed effects parameters
    betaHat <- matrix(NA, nrow = p, ncol = L)
    lambda <- rep(NA, p)

    for (r in 1:p) {
      fit_smooth <- gam(
        mum$betaTilde[r,] ~ s(argvals, bs = splines, k = (nknots + 1)),
        method = smooth_method
      )
      betaHat[r,] <- fit_smooth$fitted.values
      lambda[r] <- fit_smooth$sp # Smoothing parameter
    }

    sm <- smoothCon(
      s(argvals, bs = splines, k = (nknots + 1)),
      data = data.frame(argvals = argvals),
      absorb.cons = TRUE
    )
    S <- sm[[1]]$S[[1]] # Penalty matrix
    B <- sm[[1]]$X # Basis functions
    rm(fit_smooth, sm)
  } else {
    betaHat <- t(
      apply(
        mum$betaTilde,
        1,
        function(x) {
          gam(
            x ~ s(argvals, bs = splines, k = (nknots + 1)),
            method = smooth_method
          )$fitted.values
        }
      )
    )
  }
  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- 1:L

  # 3. Variance estimation #####################################################

  # 3.0 Early return ===========================================================

  # End the function call if no variance calculation is required
  if (var == FALSE) {
    if (!silent) {
      message(
        paste0(
          "Complete!", "\n",
          " - Use plot_fui() function to plot estimates", "\n",
          " - For more information, run the command:  ?plot_fui"
        )
      )
    }
    # Add additional return of smoothed HHat
    return(list(betaHat = betaHat, argvals = argvals, aic = mum$AIC_mat))
  }

  # At this point, the function either chooses analytic or bootstrap inference
  # Uses bootstrap in the non-analytic case

  # 3.1 Analytic inference #####################################################

  if (analytic) {
    if (!silent) print("Step 3: Inference (Analytic)")
    var_res <- var_analytic(
      mum = mum,
      betaHat = betaHat,
      data = data,
      L. = L,
      MoM = MoM,
      non_neg = non_neg,
      nknots_cov = nknots_cov,
      parallel = parallel,
      silent = silent
    )
  } else {
    if (!silent) print("Step 3: Inference (Bootstrap)")
    var_res <- var_bootstrap(mum)
  }

  return(var_res)

}
