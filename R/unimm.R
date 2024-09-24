#' Create a new "unimm" object
#'
#' The class "unimm" (and its inheritors) contain parameters for the massively
#' univariate model step. The method `fit_unimm` will be dispatched to fit the
#' mixed models.
#'
#' @param model_formula Nested character list created from a formula object in
#' lme4 formula syntax.
#' @param family GLM family of the response. Uses value fed to `fui`.
#' @param residuals Logical, indicating whether to save residuals from
#' unsmoothed LME. Uses value fed to `fui`.
#' @param caic Logical, indicating whether to calculate cAIC. Defaults to
#' \code{FALSE}.
#' @param REs Logical, indicating whether to return random effect estimates.
#' Uses value fed to `fui`.
#' @param analytic Logical, indicating whether to use the analytic inference
#' approach or bootstrap. Uses value fed to `fui`.
#'
#' @return A "unimm" object containing parameters for the univariate step.

new_unimm <- function(
  model_formula,
  family,
  residuals,
  caic,
  REs,
  analytic
) {
  unimm <- list(
    model_formula = model_formula,
    family = family,
    residuals = residuals,
    caic = caic,
    REs = REs,
    analytic = analytic
  )
  class(unimm) <- "unimm"
  return(unimm)
}

#' Create a new "unimm_conc" object
#'
#' The class "unimm" (and its inheritors) contain parameters for the massively
#' univariate model step. The "unimm_conc" class is meant for fitting concurrent
#' models. The method `fit_unimm` will be dispatched.
#'
#' @param model_formula Nested character list created from a formula object in
#' lme4 formula syntax.
#' @param family GLM family of the response. Uses value fed to `fui`.
#' @param residuals Logical, indicating whether to save residuals from
#' unsmoothed LME. Uses value fed to `fui`.
#' @param caic Logical, indicating whether to calculate cAIC. Defaults to
#' \code{FALSE}.
#' @param REs Logical, indicating whether to return random effect estimates.
#' Uses value fed to `fui`.
#' @param analytic Logical, indicating whether to use the analytic inference
#' approach or bootstrap. Uses value fed to `fui`.
#' @param func_covs Character vector of names of functional covariates.
#'
#' @return A "unimm_conc" object containing parameters for the univariate step.

new_unimm_conc <- function(
  model_formula,
  family,
  residuals,
  caic,
  REs,
  analytic,
  func_covs
) {
  unimm_conc <- new_unimm(
    model_formula,
    family,
    residuals,
    caic,
    REs,
    analytic
  )

  unimm_conc$func_covs <- func_covs
  class(unimm_conc) <- c("unimm_conc", class(unimm_conc))
  return(unimm_conc)
}
