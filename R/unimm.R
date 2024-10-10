#' Create a new "unimm" object
#'
#' The class "unimm" (and its inheritors) contain parameters for the massively
#' univariate model step. The method `fit_unimm` will be dispatched to fit the
#' mixed models.
#'
#' @param formula Formula in `lme4` formula syntax.
#' @param family Character, GLM family of the response. Passed from `fui`.
#' @param residuals Logical, indicates whether residuals are saved from
#' unsmoothed LME. Passed from `fui`.
#' @param caic Logical, indicates cAIC calculation return. Defaults to `FALSE`.
#' @param REs Logical, indicates random effect estimates return. Passed from
#' `fui`.
#' @param analytic Logical, indicates whether analytic variance will be needed.
#' Passed from `fui`.
#' @param concurrent Logical, indicates whether model fitting is concurrent.
#' @param func_covs Character vector of functional covariate names.
#'
#' @return A "unimm" object containing parameters for the univariate step. The
#' class will be "unimm_conc" if indicated by `concurrent = TRUE`.

new_unimm <- function(
  formula, family, residuals, caic, REs, analytic, concurrent
) {
  unimm <- list(
    formula = formula,
    family = family,
    residuals = residuals,
    caic = caic,
    REs = REs,
    analytic = analytic,
    concurrent = concurrent,
    func_covs = func_covs
  )

  class(unimm) <- "unimm"

  # Create a concurrent model if functional covariates are provided
  if (concurrent) {
    unimm$func_covs <- func_covs
    class(unimm) <- c("unimm_conc", class(unimm))
  }

  return(unimm)
}

#' Generic "fit_unimm"
#'
#' Dispatches class-specific methods for fitting a single point during the
#' univariate mixed model step.
#'
#' @param uni_model An object that is or inherits from the "unimm" class.
#' @param ... Additional arguments.
#'
#' @return Results depends on the dispatched method.

fit_unimm <- function(uni_model, ...) {
  UseMethod("fit_unimm")
}

#' Fit a univariate mixed model
#'
#' Fits a mixed model at location l. Part of Step 1 of the FUI approach.
#'
#' @param uni_model A `unimm` object that contains parameters for the fit.
#' @param l location to fit the model
#' @param data data frame containing all the variables in formula. Uses value
#' @param ... Additional arguments (currently ignored)
#'
#' @method fit_unimm unimm
#' @import lme4
#' @importFrom stats as.formula
#'
#' @return a list containing point estimates, variance estimates, etc.

fit_unimm.unimm <- function(uni_model, l, data, ...) {
  # Extract the data at the given point l
  form <- as.character(uni_model$formula)
  out_index <- grep(paste0("^", form[2]), names(data))
  data$Yl <- unclass(data[, out_index][, l])
  # Create a new formula
  formula <- stats::as.formula(paste0("Yl ~ ", form[3]))

  # Fit an lmer or glmer model
  fit_uni <- unimm_lmer(formula, data, uni_model$family)
  res <- unimm_outs(fit_uni, uni_model)
  res$out_index <- out_index

  return(res)
}

#' Fit a univariate mixed model for concurrent models
#'
#' Fits a mixed model at location l. Part of Step 1 of the FUI approach. Returns
#' the list of Z matrices for concurrent estimation.
#'
#' @param uni_model A `unimm` object that contains parameters for the fit.
#' @param l location to fit the model
#' @param data data frame containing all the variables in formula. Uses value
#' @param ... Additional arguments (currently ignored)
#'
#' @method fit_unimm unimm_conc
#' @import lme4
#' @importFrom stats as.formula model.matrix
#'
#' @return a list containing point estimates, variance estimates, etc. Also
#' contains `ztlist`, a list of transposed `Z` matrices.

fit_unimm.unimm_conc <- function(uni_model, l, data, ...) {
  # Extract the correct index of the functional domain
  form <- as.character(uni_model$formula)
  out_index <- grep(paste0("^", form[2]), names(data))
  data$Yl <- unclass(data[, out_index][, l])

  # Get the nonfunctional variables
  all_vars <- all.vars(uni_model$formula)
  nonfunc_vars <- all_vars[!all_vars %in% uni_model$func_vars]

  # Extract the correct indices for each functional variable
  temp <- data[nonfunc_vars]
  for (nm in uni_model$func_vars) {
    out_index <- grep(paste0("^", nm), names(dat))
    temp[[paste0(nm, "_", l)]] <- unclass(dat[, out_index][, l])
  }

  # Replace functional covariates in the formula
  formula_l <- gsub(
    paste0("(", paste0(func_vars, collapse = "|"), ")"),
    paste0("\\1_", l),
    deparse(uni_model$formula)
  )

  # Fit an lmer or glmer model
  fit_uni <- unimm_lmer(formula_l, data, uni_model$family)
  res <- unimm_outs(fit_uni, uni_model)

  # Z matrix and random effects table
  if (uni_model$analytic) {
    # AX: Add rowsums for mom == 1 case
    # colSums
    # AX: Check that the downstream casting works (e.g., matrix(x, ncol = 1))
    res$ztlist <- sapply(lme4::getME(fit_uni, "Ztlist"), function(x) colSums(x))
    varcorr_df <- as.data.frame(lme4::VarCorr(fit_uni))
    res$varcorr_df <- varcorr_df[varcorr_df$grp != "Residual", 1:3]
    res$designmat <- stats::model.matrix(fit_uni)
  }

  res$out_index <- out_index
  return(res)
}

#' Fit and return a univariate mixed model (helper)
#'
#' Helper for `fit_unimm`. Detects whether to use `lmer` or `glmer` and returns
#' the resulting model. Also helpful for producing a sample model during the
#' massively univariate step.
#'
#' @param formula Model formula.
#' @param data Data frame to use to fit the model.
#' @param family Model family. Use `lmer` if "gaussian", `glmer` otherwise.
#'
#' @import lme4
#'
#' @return an `lme4` model

unimm_lmer <- function(formula, data, family) {
  if (family == "gaussian") {
    suppressMessages(
      lme4::lmer(
        formula = formula,
        data = data,
        control = lme4::lmerControl(
          optimizer = "bobyqa", optCtrl = list(maxfun = 5000)
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
          optimizer = "bobyqa", optCtrl = list(maxfun = 5000)
        )
      )
    )
  }
}

#' Get relevant features of univariate models
#'
#' Helper for `fit_unimm` that returns various qualities of the univariate fit
#' produced by `unimm_lmer` (also a helper).
#'
#' @param fit_uni An `lme4` object corresponding to the fit at some point on the
#' functional domain.
#' @param uni_model A "unimm" class object with parameters of the fit
#'
#' @return A list of relevant features of `fit_uni`
#'
#' @import lme4
#' @importFrom cAIC4 cAIC
#' @importFrom stats residuals AIC BIC

unimm_outs <- function(fit_uni, uni_model) {
  # Fixed effects estimates
  betaTilde <- lme4::fixef(fit_uni)

  # Initialize returned parameters
  randeffs <- aic_met <- resids <- NA

  # these are residuals from including the random effects (i.e., with BLUPs),
  # not JUST from fixed effects
  # Compare with nlme::lme(); 2 columns of residuals in lme: mod$residuals
  if(uni_model$residuals) resids <- as.numeric(residuals(fit_uni))
  if(uni_model$caic) aic_met <- as.numeric(cAIC4::cAIC(fit_uni)$caic)[1]
  # random effects
  if(uni_model$REs) randeffs <- lme4::ranef(fit_uni)
  varcorr <- as.data.frame(lme4::VarCorr(fit_uni))

  # Setup returned results
  res <- list(
    betaTilde = betaTilde,
    group = varcorr[1, 1],
    aic = stats::AIC(fit_uni),
    bic = stats::BIC(fit_uni),
    residuals = resids,
    caic = aic_met,
    randeffs = randeffs
  )

  if (!uni_model$analytic)
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
