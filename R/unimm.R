#' Generic "unimm" model fitting
#'
#' Dispatches class-specific methods for fitting a single point during the
#' univariate mixed model step.
#'
#' @param fmm An object that is or inherits from the "fastFMM" class.
#' @param ... Additional arguments.
#'
#' @return Results depends on the dispatched method.
#' @export

unimm <- function(fmm, ...) {
  UseMethod("unimm")
}

#' Fit a univariate mixed model
#'
#' Fits a mixed model at location l. Part of Step 1 of the FUI approach.
#'
#' @param fmm A `fastFMM` object that contains parameters for the fit.
#' @param l location to fit the model
#'
#' @method unimm fastFMM
#' @import lme4
#' @importFrom stats as.formula
#'
#' @return a list containing point estimates, variance estimates, etc.
#' @export

unimm.fastFMM <- function(fmm, l) {
  # Extract the data at the given point l
  form <- as.character(fmm$formula)
  out_index <- fmm$out_index
  data <- fmm$data
  data$Yl <- unclass(data[, out_index][, l])
  # Create a new formula
  formula <- stats::as.formula(paste0("Yl ~ ", form[3]))

  # Fit an lmer or glmer model
  fit_uni <- unimm_lmer(formula, data, fmm$family)
  res <- unimm_outs(fit_uni, fmm)

  return(res)
}

#' Fit a univariate mixed model for concurrent models
#'
#' Fits a mixed model at location l. Part of Step 1 of the FUI approach. Returns
#' the list of Z matrices for concurrent estimation.
#'
#' @param fmm A `fastFMM` object that contains parameters for the fit.
#' @param l location to fit the model
#' @param MoM indicator of type of MoM estimator. Coerced to `MoM == 1` to
#' prevent storage of all Z matrices.
#'
#' @method unimm fastFMMconc
#' @return a list containing point estimates, variance estimates, etc. Also
#' contains `ztlist`, a list of transposed `Z` matrices.
#'
#' @import lme4
#' @importFrom stats as.formula model.matrix
#' @export

unimm.fastFMMconc <- function(fmm, l) {
  # Extract the correct index of the functional domain
  form <- as.character(fmm$formula)
  data <- fmm$data
  out_index <- fmm$out_index
  temp <- data[, -out_index]

  # Get the nonfunctional variables
  all_vars <- all.vars(fmm$formula)
  fun_covariates <- fmm$fun_covariates
  nonfun_covariates <- all_vars[
    !all_vars %in% c(fmm$fun_covariates, form[2])
  ]
  # Extract the correct indices for each functional variable
  temp <- temp[nonfun_covariates]
  temp$Yl <- unclass(data[, out_index][, l])
  for (nm in fmm$fun_covariates) {
    fun_index <- grep(paste0("^", nm), names(data))
    temp[[paste0(nm, "_", l)]] <- unclass(data[, fun_index][, l])
  }

  # Replace functional covariates in the formula
  formula_l <- paste0(
    "Yl ~ ",
    gsub(
      # Regex of the functional covariate names
      paste0("(", paste0(fun_covariates, collapse = "|"), ")"),
      paste0("\\1_", l),
      form[3]
    )
  )
  formula_l <- as.formula(formula_l)

  # Fit an lmer or glmer model
  fit_uni <- unimm_lmer(formula_l, temp, fmm$family)
  res <- unimm_outs(fit_uni, fmm)

  # Z matrix and random effects table
  if (fmm$analytic) {
    res$ztlist <- sapply(lme4::getME(fit_uni, "Ztlist"), function(x) colSums(x))
    varcorr_df <- as.data.frame(lme4::VarCorr(fit_uni))
    res$varcorr_df <- varcorr_df[varcorr_df$grp != "Residual", 1:3]
    res$designmat <- stats::model.matrix(fit_uni)
  }

  return(res)
}

#' Fit and return a univariate mixed model (helper)
#'
#' Helper for `unimm`. Detects whether to use `lmer` or `glmer` and returns
#' the resulting model. Also helpful for producing a sample model during the
#' massively univariate step.
#'
#' @param formula Model formula.
#' @param data Data frame of observations to fit
#' @param family Character, family of model
#'
#' @import lme4
#'
#' @return an `lme4` model chosen between `lmer` or `glmer`, as appropriate.
#' @export

unimm_lmer <- function(formula, data, family) {
  if (family == "gaussian") {
    suppressMessages(
      lme4::lmer(
        formula = formula,
        data = data,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 5000),
          check.rankX = "ignore"
        )
      )
    )
  } else {
    suppressMessages(
      lme4::glmer(
        formula = formula,
        data = data,
        family = family,
        control = lme4::glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 5000),
          check.rankX = "ignore"
        )
      )
    )
  }
}

#' Get relevant features of univariate models
#'
#' Helper for `unimm` that returns various qualities of the univariate fit
#' produced by `unimm_lmer` (also a helper).
#'
#' @param fit_uni An `lme4` object corresponding to the fit at some point on the
#' functional domain.
#' @param fmm A "fastFMM" class object with parameters of the fit
#'
#' @return A list of relevant features of `fit_uni`
#'
#' @import lme4
#' @importFrom cAIC4 cAIC
#' @importFrom stats residuals AIC BIC
#' @export

unimm_outs <- function(fit_uni, fmm) {
  # Fixed effects estimates
  betaTilde <- lme4::fixef(fit_uni)

  # Initialize returned parameters
  randeffs <- aic_met <- residuals <- NA

  # these are residuals from including the random effects (i.e., with BLUPs),
  # not JUST from fixed effects
  # Compare with nlme::lme(); 2 columns of residuals in lme: mod$residuals
  if(fmm$residuals) residuals <- as.numeric(stats::residuals(fit_uni))
  if(fmm$caic) aic_met <- as.numeric(cAIC4::cAIC(fit_uni)$caic)[1]
  # random effects
  if(fmm$randeffs) randeffs <- lme4::ranef(fit_uni)
  varcorr <- as.data.frame(lme4::VarCorr(fit_uni))

  # Setup returned results
  res <- list(
    betaTilde = betaTilde,
    group = varcorr[1, 1],
    aic = stats::AIC(fit_uni),
    bic = stats::BIC(fit_uni),
    residuals = residuals,
    caic = aic_met,
    randeffs = randeffs
  )

  if (!fmm$analytic)
    return(res)

  # Need additional info for analytic variance calculation
  # Extract variance/covariance estimates
  var_random <- varcorr[, 4]
  # Variance of random components
  ind_var <- which(is.na(varcorr[, 3]) & varcorr[, 1] != "Residual")
  names(var_random)[ind_var] <- paste0(
    "var.", varcorr[ind_var, 1], ".", varcorr[ind_var,2]
  )

  # variance of the residual components
  names(var_random)[which(varcorr[, 1] == "Residual")] <- "var.Residual"

  # covariance of random components
  names(var_random)[which(!is.na(varcorr[,3]))] <- paste0(
    "cov.",
    varcorr$grp[which(!is.na(varcorr[, 3]))], ".",
    varcorr$var1[which(!is.na(varcorr[,3]))], ".",
    varcorr$var2[which(!is.na(varcorr[,3]))]
  )

  # SE of fixed effects
  res$var_random <- var_random
  res$se_mat <- summary(fit_uni)$coefficients[, 2]

  return(res)
}
