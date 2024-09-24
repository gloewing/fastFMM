#' Generic "massmm" method
#'
#' Create a new "massmm" object based on the class of the univariate model.
#'
#' @param obj An object that is or inherits from the "unimm" class.
#' @param ... Additional arguments.
#'
#' @return Results depends on the dispatched method.

massmm <- function(obj, ...) {
  UseMethod("massmm")
}

#' Create a new "massmm" object
#'
#' @param obj A "unimm" object that provides parameters for the fit.
#' @param res A list of objects returned from `fit_unimm` after the massively
#' univariate step
#' @param ... Additional arguments (ignored)
#'
#' @return A "massmm" object.

massmm.unimm <- function(res, unimm_obj) {
  # AX: Consider initializing the returned list here
  # e.g., massmm <- list(...)

  # Obtain betaTilde, fixed effects estimates
  betaTilde <- t(do.call(rbind, lapply(res, '[[', 1)))
  colnames(betaTilde) <- argvals

  # Obtain residuals, AIC, BIC, and random effects estimates (analytic)
  ## AIC/BIC
  mod_aic <- do.call(c, lapply(res, '[[', 'aic'))
  mod_bic <- do.call(c, lapply(res, '[[', 'bic'))
  mod_caic <- do.call(c, lapply(res, '[[', 'caic'))
  AIC_mat <- cbind(mod_aic, mod_bic, mod_caic)
  colnames(AIC_mat) <- c("AIC", "BIC", "cAIC")

  ## Store residuals if resids == TRUE
  resids <- NA
  if (unimm_obj$residuals)
    resids <- suppressMessages(
      lapply(res, '[[', 'residuals') %>% dplyr::bind_cols()
    )

  # If not analytic, stop here
  if (!unimm_obj$analytic) {
    massmm <- list(
      betaTilde = betaTilde,
      AIC_mat = AIC_mat,
      resids = resids,
      analytic = unimm_obj$analytic
    )
    class(massmm) <- "massmm"
    return(massmm)
  }

  ## random effects
  # Initialize random effects and SE matrix
  # AX: Check if this line can be removed entirely
  randEff <- ifelse(
    unimm_obj$REs,
    suppressMessages(
      simplify2array(
        lapply(lapply(res, '[[', 're_df'), function(x) as.matrix(x[[1]]))
      )
    ),
    NULL
  )
  se_mat <- suppressMessages(do.call(cbind, lapply(res, '[[', 9)))

  var_random <- t(do.call(rbind, lapply(massmm, '[[', 'var_random')))
  sigmaesqHat <- var_random["var.Residual", , drop = FALSE]
  sigmausqHat <- var_random[
    which(rownames(var_random) != "var.Residual"), , drop = FALSE
  ]

  # Fit a fake model to obtain a design matrix
  # Fit a fake model to obtain design matrix
  data$Yl <- unclass(data[,out_index][, 1])
  if (family == "gaussian") {
    fit_uni <- suppressMessages(
      lmer(
        formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
        data = data,
        control = lmerControl(optimizer = "bobyqa")
      )
    )
  } else {
    fit_uni <- suppressMessages(
      glmer(
        formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
        data = data,
        family = family,
        control = glmerControl(optimizer = "bobyqa")
      )
    )
  }

  # Design matrix
  designmat <- stats::model.matrix(fit_uni)
  name_random <- as.data.frame(VarCorr(fit_uni))[
    which(!is.na(as.data.frame(VarCorr(fit_uni))[, 3])), 3
  ]

  # Names of random effects
  RE_table <- as.data.frame(VarCorr(fit_uni))
  ranEf_grp <- RE_table[, 1]
  RE_table <- RE_table[RE_table$grp != "Residual", 1:3]
  ranEf_grp <- ranEf_grp[ranEf_grp != "Residual"]
  ztlist <- sapply(getME(fit_uni, "Ztlist"), t)

  # Check if group contains ":" (hierarchical structure) that requires the
  # group to be specified

  group <- massmm[[1]]$group ## group name in the data
  if (grepl(":", group, fixed = TRUE)) {
    if (is.null(subj_ID)) {
      # Assumes the ID name is to the right of the ":"
      group <- str_remove(group, ".*:")
    } else {
      # Use user-specified ID if it exists
      group <- subj_ID
    }
  }

  if (is.null(subj_ID))
    subj_ID <- group

  # Condition for using G_estimate_randint
  randint_flag <- I(
    length(fit_uni@cnms) == 1 &
      length(fit_uni@cnms[[group]]) == 1 &
      fit_uni@cnms[[group]][1] == "(Intercept)"
  )

  # Check for what needs to be returned
  # This will get passed to smoothing
  massmm <- list(
    betaTilde = betaTilde,
    AIC_mat = AIC_mat,
    RE_table = RE_table,
    ...,
    randint_flat = randint_flag
  )
  class(massmm) <- "massmm"
  return(massmm)
}

massmm.unimm_conc <- function(obj, res, ...) {

}
