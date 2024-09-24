#' Generic "fit_unimm"
#'
#' Dispatches class-specific methods for fitting a single point during the
#' univariate mixed model step.
#'
#' @param obj An object that is or inherits from the "unimm" class.
#' @param ... Additional arguments.
#'
#' @return Results depends on the dispatched method.

fit_unimm <- function(obj, ...) {
  UseMethod("fit_unimm")
}

#' Fit a univariate mixed model
#'
#' Fits a mixed model at location l. Part of Step 1 of the FUI approach.
#'
#' @param obj A `unimm` object that contains parameters for the fit.
#' @param l location to fit the model
#' @param data data frame containing all the variables in formula. Uses value
#' @param ... Additional arguments (currently ignored)
#'
#' @method fit_unimm unimm
#' @import lme4
#' @importFrom stats as.formula AIC BIC
#' @importFrom cAIC4 cAIC
#'
#' @return a list containing point estimates, variance estimates, etc.

fit_unimm.unimm <- function(obj, l, data, ...) {
  # Extract the data at the given point l
  form <- obj$model_formula
  out_index <- grep(paste0("^", form[2]), names(data))
  data$Yl <- unclass(data[, out_index][, l])

  # Fit the model as a GLM if the family is not Gaussian
  if(obj$family == "gaussian"){
    fit_uni <- suppressMessages(
      lme4::lmer(
        formula = stats::as.formula(paste0("Yl ~ ", form[3])),
        data = data,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 5000)
        )
      )
    )
  } else {
    fit_uni <- suppressMessages(
      lme4::glmer(
        formula = stats::as.formula(paste0("Yl ~ ", form[3])),
        data = data,
        family = family,
        control = lme4::glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 5000)
        )
      )
    )
  }

  # Fixed effects estimates
  betaTilde <- lme4::fixef(fit_uni)

  # Initialize returned parameters
  re_df <- aic_met <- resids <- NA

  # these are residuals from including the random effects (i.e., with BLUPs),
  # not JUST from fixed effects
  # Compare with nlme::lme(); 2 columns of residuals in lme: mod$residuals
  if(obj$residuals) resids <- as.numeric(residuals(fit_uni))
  if(obj$caic) aic_met <- as.numeric(cAIC4::cAIC(fit_uni)$caic)[1]
  # random effects
  if(obj$REs) re_df <- lme4::ranef(fit_uni)
  varcorr <- as.data.frame(lme4::VarCorr(fit_uni))

  # Setup returned results
  unimm_res <- list(
    betaTilde = betaTilde,
    group = varcorr[1, 1],
    aic = stats::AIC(fit_uni),
    bic = stats::BIC(fit_uni),
    residuals = resids,
    caic = aic_met,
    re_df = re_df
  )

  # Return now if not analytic
  if (!obj$analytic)
    return(unimm_res)

  # Continue with calculations of SE if analytic

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
  names(var_random)[
    which(!is.na(as.data.frame(lme4::VarCorr(fit_uni))[,3]))
  ] <- paste0(
    "cov.",
    varcorr$grp[
      which(!is.na(as.data.frame(lme4::VarCorr(fit_uni))[,3]))
    ],
    ".",
    varcorr$var1[
      which(!is.na(as.data.frame(lme4::VarCorr(fit_uni))[,3]))
    ],
    ".",
    varcorr$var2[
      which(!is.na(as.data.frame(lme4::VarCorr(fit_uni))[,3]))
    ]
  )

  # SE of fixed effects
  unimm_res$var_random <- var_random
  unimm_res$se_mat <- summary(fit_uni)$coefficients[, 2]

  # Analytic results
  unimm_res
}

# Verify the arguments
fit_unimm.unimm_conc <- function(obj, l, data, ...) {
  # Extract the correct index of the functional domain
  out_index <- grep(paste0("^", obj$model_formula[2]), names(data))
  data$Yl <- unclass(data[,out_index][,l])
  out_index <- grep(paste0("^", model_formula[2]), names(data))
  cov_index <- grep(paste0("^", fun_cov), names(data))

  data$Yl <- unclass(data[, out_index][, l])
  data$fun_cov <- unclass(data[, cov_index][, l])

  if (family == "gaussian") {
    fit_uni <- suppressMessages(
      lmer(
        formula = stats::as.formula(paste0("Yl ~ fun_cov + ", model_formula[3])),
        data = data,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 5000)
        )
      )
    )
  } else {
  fit_uni <- suppressMessages(
      glmer(
        formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
        data = data,
        family = family,
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5000))
      )
    )
  }

  # Fixed effects estimates
  betaTilde <- lme4::fixef(fit_uni)
  re_df <- aic_met <- resids <- NA

  # these are residuals from including the random effects (i.e., with BLUPs),
  # not JUST from fixed effects
  # Compare with nlme::lme(); 2 columns of residuals in lme: mod$residuals
  if(residuals) resids <- as.numeric(residuals(fit_uni))
  if(caic) aic_met <- as.numeric(cAIC4::cAIC(fit_uni)$caic)[1]
  # random effects
  if(REs) re_df <- ranef(fit_uni)

  if(analytic == TRUE) {
    varcorr <- as.data.frame(VarCorr(fit_uni))
    # Extract variance/covariance estimates
    var_random <- varcorr[,4]

    # Variance of random components
    ind_var <- which(is.na(varcorr[,3]) & varcorr[,1] != "Residual")
    names(var_random)[ind_var] <- paste0("var.",varcorr[ind_var,1],".",varcorr[ind_var,2])

    # variance of the residual components
    names(var_random)[which(varcorr[,1] == "Residual")] <- "var.Residual"

    # covariance of random components
    names(var_random)[which(!is.na(as.data.frame(VarCorr(fit_uni))[,3]))] <-
      paste0("cov.",
             varcorr$grp[which(!is.na(as.data.frame(VarCorr(fit_uni))[,3]))], ".",
             varcorr$var1[which(!is.na(as.data.frame(VarCorr(fit_uni))[,3]))], ".",
             varcorr$var2[which(!is.na(as.data.frame(VarCorr(fit_uni))[,3]))])
    se_mat <- summary(fit_uni)$coefficients[,2] ## se of fixed effects

    # Z matrix
    ztlist <- sapply(getME(fit_uni, "Ztlist"), function(x) t(x) )

    return(
      list(
        betaTilde = betaTilde,
        group = varcorr[1,1],
        aic = stats::AIC(fit_uni),
        bic = stats::BIC(fit_uni),
        residuals = resids,
        caic = aic_met,
        re_df = re_df,
        var_random = var_random,
        se_mat = se_mat,
        ztlist = ztlist
      )
    )
  } else {
    return(
      list(
        betaTilde = betaTilde,
        group = as.data.frame(VarCorr(fit_uni))[1,1],
        aic = stats::AIC(fit_uni),
        bic = stats::BIC(fit_uni),
        residuals = resids,
        caic = aic_met,
        re_df = re_df,
        ztlist = ztlist
      )
    )
  }
}
