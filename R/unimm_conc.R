#' Fit a univariate mixed model
#'
#' Fits a mixed model at location l. Part of Step 1 of the FUI approach.
#' Unlike the `unimm` function for non-concurrent FMM, this function returns
#' the Z matrix at location l.
#'
#' @param l location to fit the model
#' @param data data frame containing all the variables in formula. Uses value
#' fed to `fui`.
#' @param model_formula Character version of a two-sided formula object in lme4
#' formula syntax, produced within `fui`.
#' @param family GLM family of the response. Uses value fed to `fui`.
#' @param residuals Logical, indicating whether to save residuals from
#' unsmoothed LME. Uses value fed to `fui`.
#' @param caic Logical, indicating whether to calculate cAIC. Defaults to \code{FALSE}.
#' @param REs Logical, indicating whether to return random effect estimates.
#' Uses value fed to `fui`.
#' @param analytic Logical, indicating whether to use the analytic inference
#' approach or bootstrap. Uses value fed to `fui`.
#' @param cov Character, name of the functional covariate
#'
#' @import lme4
#'
#' @return a list containing point estimates, variance estimates, etc.

# TODO: generalize to multiple functional cov

unimm_conc <- function(
  l,
  data,
  model_formula,
  family,
  residuals,
  caic,
  REs,
  analytic,
  fun_cov
) {

  out_index <- grep(paste0("^", model_formula[2]), names(data))
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
