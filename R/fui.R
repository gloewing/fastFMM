#' Fast Univariate Inference for Longitudinal Functional Models
#'
#' Fit a function-on-scalar regression model for longitudinal
#' functional outcomes and scalar predictors using the Fast Univariate
#' Inference (FUI) approach (Cui et al. 2022; Loewinger et al. 2024).
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
#' by Cui et al. (2022). For more information on the method of moments estimator 
#' applied in step 3, see Loewinger et al. (2024).
#'
#' @param formula Two-sided formula object in lme4 formula syntax.
#' The difference is that the response need to be specified as a matrix
#' instead of a vector. Each column of the matrix represents one location
#' of the longitudinal functional observations on the domain.
#' @param data A data frame containing all variables in formula
#' @param family GLM family of the response. Defaults to \code{gaussian}.
#' @param var Logical, indicating whether to calculate and return variance
#' of the coefficient estimates. Defaults to `TRUE`.
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
#' @param n_cores Number of cores for parallelization. Defaults to 1.
#' @param caic Logical, indicating whether to calculate cAIC. Defaults to
#' \code{FALSE}.
#' @param randeffs Logical, indicating whether to return random effect estimates.
#' Defaults to \code{FALSE}.
#' @param non_neg 0 - no non-negativity constrains, 1 - non-negativity
#' constraints on every coefficient for variance, 2 - non-negativity on
#' average of coefficents for 1 variance term. Defaults to 0.
#' @param MoM Method of moments estimator. Defaults to 1.
#' @param concurrent Logical, indicates whether to fit a concurrent model.
#' Defaults to \code{FALSE}.
#' @param impute_outcome Logical, indicates whether to impute missing outcome
#' values with FPCA. Defaults to \code{FALSE}. Use with caution as the
#' downstream effects are not tested.
#' @param override_zero_var Logical, indicates whether to proceed with model
#' fitting if columns have zero variance. Suggested for cases where individual
#' columns have zero variance but interactions have non-zero variance. Defaults
#' to `FALSE`.
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
#' @references Loewinger, G., Cui, E., Lovinger, D., Pereira, F. (2024). A 
#' Statistical Framework for Analysis of Trial-Level Temporal Dynamics in 
#' Fiber Photometry Experiments. \emph{eLife}, 95802.
#'
#' @export
#'
#' @import lme4
#' @import parallel
#' @import magrittr
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
  n_cores = 1,
  caic = FALSE,
  randeffs = FALSE,
  non_neg = 0,
  MoM = 1,
  concurrent = FALSE,
  impute_outcome = FALSE,
  override_zero_var = FALSE
) {

  # 0. Setup ###################################################################

  # 0.0 Argument consistency checks ============================================

  # If doing parallel computing, set up the number of cores
  if (parallel & !is.integer(n_cores)) {
    n_cores <- as.integer(round(parallel::detectCores() * 0.75))
    if (n_cores < 2)
      warning("Only 1 core detected for parallelization.")
    if (!silent)
      message("Cores used for parallelization: ", n_cores)
  }

  # For non-Gaussian family, manually set variance to bootstrap inference
  if (family != "gaussian") {
    if (analytic & !silent) { # Notify user of conflict
      message(
        "Analytic variance is not supported for non-Gaussian models. ",
        "Variance calculation will be done through bootstrap."
      )
    }
    analytic <- FALSE
  }

  # 0.1 Concurrent model checks ================================================

  # Detect functional covariates by matching length to L
  model_formula <- as.character(formula)
  out_index <- grep(paste0("^", model_formula[2]), names(data))
  L <- length(out_index)
  x_names <- all.vars(formula)[-1]
  x_classes <- sapply(x_names, function(x) class(data[[x]]))
  x_ncols <- sapply(
    x_names,
    function(x) length(grep(paste0("^", x), names(data)))
  )
  fun_covariates <- NULL # redundancy
  fun_covariates <- unique(c(x_names[x_classes == "AsIs"], x_names[x_ncols == L]))
  fun_exists <- length(fun_covariates) > 0
  if (fun_exists)
    message("Functional covariate(s): ", paste0(fun_covariates, collapse = ", "))
  # Check for inconsistencies with user-set concurrence argument
  if (concurrent & !fun_exists) {
    stop(
      "No functional covariates found for concurrent model fitting.", "\n",
      "Incoporate functional covariates or set `concurrent = F` in args."
    )
  } else if (!concurrent & fun_exists) {
    warning(
      "Functional covariates detected: ",
      paste0(fun_covariates, collapse = ", "), "\n",
      "Execution will continue, but consider setting `concurrent = TRUE`."
    )
  }

  # Check for the MoM estimator and coerce to 1
  if (concurrent & MoM == 2) {
    warning(
      "MoM = 2 is currently not supported for concurrent models. ",
      "Calculation will proceed with MoM = 1."
    )
    MoM <- 1
  }

  # 0.0 Identifiability checks ==================================================

  all_vars <- all.vars(formula)
  out_index <- grep(
    paste0("^(", paste0(all_vars, collapse = "|"), ")"),
    names(data)
  )
  temp <- data[, out_index]
  # Coerce characters, like IDs, to numerics
  temp <- data.frame(
    lapply(
      temp,
      function(col) {
        if (is.numeric(col))
          return(col)
        as.numeric(as.factor(col))
      }
    )
  )
  col_var <- Rfast::colVars(as.matrix(temp))
  col_var_zero <- which(col_var == 0)

  if (length(col_var_zero) > 0) {
    msg <- paste0(
      "Columns with zero variance: ",
      paste0(names(temp)[col_var_zero], collapse = ", "), "\n",
      "Model-fitting cannot continue due to non-identifiability."
    )
    ifelse(override_zero_var, warning(msg), stop(msg))
  }

  # 0.2 Create the reference object ============================================

  fmm_params <- list(
    formula = formula,
    data = data,
    subj_id = subj_id,
    argvals = argvals,
    family = family,
    residuals = residuals,
    caic = caic,
    randeffs = randeffs,
    var = var,
    analytic = analytic
  )

  if (!concurrent) {
    fmm <- do.call(new_fastFMM, fmm_params)
  } else {
    fmm_params$fun_covariates <- fun_covariates
    fmm <- do.call(new_fastFMMconc, fmm_params)
  }

  # 0.3 Impute missing values ==================================================

  out_index <- fmm$out_index
  missing_rows <- which(rowSums(is.na(data[, out_index])) != 0)

  # Fill in missing values of functional outcome using FPCA
  # rows with missing outcome values
  missing_rows <- which(rowSums(is.na(data[, out_index])) != 0 )

  if (length(missing_rows) != 0) {
    if(analytic & impute_outcome) {
      message(
        paste(
          "Imputing", sum(is.na(data[, out_index])),
          "values in functional response with longitudinal functional PCA"
        )
      )
      if (length(out_index) != 1) {
        nknots_fpca <- min(round(length(out_index) / 2), 35)
        if (is.null(argvals) | analytic)
          argvals <- 1:length(out_index)
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
          "Removing", length(missing_rows),
          "rows with missing functional outcome values.", "\n",
          "To impute missing outcome values with FPCA, set fui() argument: \n",
          "impute_outcome = TRUE"
        )
      )
      # remove data with missing rows
      data <- data[-missing_rows, ]
    }
    fmm$data <- data
  }

  # 1. Massively univariate mixed models #######################################

  if (!silent)
    print("Step 1: Fit Massively Univariate Mixed Models")

  # Create a list of univariate models ("massively univariate")
  mum <- massmm(fmm, parallel, n_cores)

  # 2. Smoothing ###############################################################

  if (!silent) print("Step 2: Smoothing")

  # Penalized splines smoothing and extract components (analytic)
  # Number of knots for regression coefficients
  nknots <- min(round(L / 2), nknots_min)
  # Number of knots for covariance matrix
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- min(round(L / 2), 35)

  # Reset argvals
  argvals <- fmm$argvals

  # Smoothing parameter, spline basis, penalty matrix (analytic)
  # Setup variables
  lambda <- S <- B <- NULL

  # Smooth coefficient estimates
  HHat <- t(
    apply(
      mum$sigmausqHat, 1,
      function(b) stats::smooth.spline(x = argvals, y = b)$y
    )
  )
  ind_var <- which(grepl("var", rownames(HHat)) == TRUE)
  HHat[ind_var, ][which(HHat[ind_var, ] < 0)] <- 0

  betaTilde <- mum$betaTilde

  if (analytic) {
    p <- nrow(betaTilde) # Number of fixed effects parameters
    betaHat <- matrix(NA, nrow = p, ncol = L)
    lambda <- rep(NA, p)

    # NB: although s() is loaded from mgcv, mgcv::s will break.
    # Spacedman describes this here: stackoverflow.com/a/20694106
    # Solution: Import mgcv in full and be careful of collisions in future
    for (r in 1:p) {
      fit_smooth <- mgcv::gam(
        betaTilde[r,] ~ s(argvals, bs = splines, k = nknots + 1),
        method = smooth_method
      )
      betaHat[r,] <- fit_smooth$fitted.values
      lambda[r] <- fit_smooth$sp # Smoothing parameter
    }

    sm <- mgcv::smoothCon(
      s(argvals, bs = splines, k = nknots + 1),
      data = data.frame(argvals = argvals),
      absorb.cons = TRUE
    )
    S <- sm[[1]]$S[[1]] # Penalty matrix
    B <- sm[[1]]$X # Basis functions
    rm(fit_smooth, sm)
  } else {
    betaHat <- t(
      apply(
        betaTilde,
        1,
        function(x) {
          mgcv::gam(
            x ~ s(argvals, bs = splines, k = nknots + 1),
            method = smooth_method
          )$fitted.values
        }
      )
    )
  }
  rownames(betaHat) <- rownames(betaTilde)
  rm(betaTilde)
  colnames(betaHat) <- 1:L

  # Save a convenient list to pass to variance calculation
  smoothed <- list(
    betaHat = betaHat,
    HHat = HHat,
    S = S,
    B = B,
    lambda = lambda
  )

  # 3. Variance estimation #####################################################

  # 3.0 Early return ===========================================================

  # End the function call if no variance calculation is required
  if (!var) {
    if (!silent) {
      message(
        paste0(
          "Complete!", "\n",
          " - Use plot_fui() function to plot estimates", "\n",
          " - For more information, run the command:  ?plot_fui"
        )
      )
    }
    # AX: Add additional return of smoothed HHat

    return(
      list(
        betaHat = smoothed$betaHat,
        HHat = smoothed$HHat,
        argvals = argvals,
        aic = mum$AIC_mat
        )
      )
  }

  # At this point, the function either chooses analytic or bootstrap inference
  # Uses bootstrap in the non-analytic case

  # 3.1 Analytic inference #####################################################

  if (analytic) {
    if (!silent) print("Step 3: Inference (Analytic)")
    var_res <- var_analytic(
      fmm,
      mum,
      smoothed,
      MoM,
      non_neg,
      nknots_cov,
      seed,
      parallel,
      n_cores,
      silent
    )
  } else {
    if (!silent) print("Step 3: Inference (Bootstrap)")
    var_res <- var_bootstrap(mum)
  }

  # AX: Can't really remove this at any other point
  if (!design_mat) var_res$designmat <- NULL


  return(var_res)
}
