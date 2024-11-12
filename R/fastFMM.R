#' Create a new "fastFMM" object
#'
#' The class "fastFMM" (and its inheritors) contain parameters for fast
#' univariate inference. The basic "fastFMM" object is equipped for
#' non-concurrent model fitting.
#'
#' Object creation populates fields relevant to later steps, such as the
#' location of the functional domain in the data frame.
#'
#' @param formula Formula in `lme4` formula syntax.
#' @param data Data frame to fit.
#' @param subj_id Character, name of variable containing IDs. Paased from `fui`
#' @param argvals List of points to fit on the functional domain. Only applies
#' for the bootstrap case (i.e., `analytic = FALSE`).
#' @param family Character, GLM family of the response. Passed from `fui`.
#' @param residuals Logical, indicates whether residuals are saved from
#' unsmoothed LME. Passed from `fui`.
#' @param caic Logical, indicates cAIC calculation return. Defaults to `FALSE`.
#' @param randeffs Logical, indicates whether random effect estimates are returned.
#' Passed from `fui`.
#' @param var Logical, indicates whether to calculate variance. Passed from
#' `fui`.
#' @param analytic Logical, indicates whether variance will be calculated
#' analytically. Passed from `fui`.
#'
#' @return A "fastFMM" object containing parameters for fitting a functional
#' mixed model using the FUI scheme. The object contains each of the passed
#' args (`formula, data, ..., analytic`), with the exception of `var`.
#' Additional entries returned are:
#' \enumerate{
#' \item `out_index`: locations in `data` where functional domain exists
#' \item `argvals`: either `argvals` as passed (if bootstrap) or a vector `1:L`
#' where `L` is the size of the functional domain
#' }
#' @export

new_fastFMM <- function(
  formula, data, subj_id, argvals, family, residuals, caic, randeffs, var, analytic
) {
  # Create basic fields of the model object
  fmm <- list(
    formula = formula,
    data = data,
    subj_id = subj_id,
    family = family,
    residuals = residuals,
    caic = caic,
    randeffs = randeffs,
    analytic = analytic
  )

  ### Populate other parameters
  # Check validity of formula
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)

  # Stop if there are column names with "." to avoid issues with G, H
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

  # Set out_index and dimension of functional domain (L)
  fmm$out_index <- out_index <- grep(paste0("^", model_formula[2]), names(data))
  # out_index may refer to a matrix column or multiple columns
  if (length(out_index) != 1) {
    # Multiple columns
    L <- length(out_index)
  } else {
    # Matrix column
    L <- ncol(data[, out_index])
  }

  # Check argvals and set appropriately
  if (analytic & !is.null(argvals) & var)
    message(
      "'argvals' argument is currently only supported for bootstrap.", "\n",
      "Overwriting argvals to fit model to ALL points on functional domain"
    )

  if (is.null(argvals) | analytic) {
    fmm$argvals <- 1:L
  } else {
    if (max(argvals) > L)
      stop(
        "Maximum index specified in argvals is greater than ",
        "the size of the functional domain for the outcome."
      )
    fmm$argvals <- argvals
    # AX: May need to set L here for boostrap
    # L <- length(argvals)
  }

  # Set the class and return the object
  class(fmm) <- "fastFMM"
  return(fmm)
}

#' Create a new "fastFMMconc" object
#'
#' Create an object that contains parameters for fast univariate inference for
#' concurrent models.
#'
#' @param formula Formula in `lme4` formula syntax.
#' @param data Data frame to fit.
#' @param argvals List of points to fit on the functional domain. Only applies
#' for the bootstrap case (i.e., `analytic = FALSE`).
#' @param family Character, GLM family of the response. Passed from `fui`.
#' @param residuals Logical, indicates whether residuals are saved from
#' unsmoothed LME. Passed from `fui`.
#' @param caic Logical, indicates cAIC calculation return. Defaults to `FALSE`.
#' @param randeffs Logical, indicates whether random effect estimates are returned.
#' Passed from `fui`.
#' @param var Logical, indicates whether to calculate variance. Passed from
#' `fui`.
#' @param analytic Logical, indicates whether variance will be calculated
#' analytically. Passed from `fui`.
#' @param fun_covariates Character vector of functional covariate names.
#'
#' @return A "fastFMMconc" object containing parameters to fit a concurrent
#' mixed model using the FUI scheme. This function is called within `fui` if
#' indicated by `concurrent = TRUE`. Fields are shared with `fastFMM` objects,
#' with the addition of `fun_covariates`.
#' @export

new_fastFMMconc <- function(
  formula,
  data,
  subj_id,
  argvals,
  family,
  residuals,
  caic,
  randeffs,
  var,
  analytic,
  fun_covariates
) {
  # Creates a basic non-concurrent model object
  fmm <- new_fastFMM(
    formula,
    data,
    subj_id,
    argvals,
    family,
    residuals,
    caic,
    randeffs,
    var,
    analytic
  )
  # Add field for functional covariates
  fmm$fun_covariates <- fun_covariates
  class(fmm) <- c("fastFMMconc", "fastFMM")
  return(fmm)
}

#' Printing the fastFMM object
#'
#' Simple method to look at the relevant fields.
#'
#' @param fmm Object to print
#'
#' @return Basic list of entries
#' @method print fastFMM
#' @export

print.fastFMM <- function(fmm) {
  paste0(
    "Formula: ", fmm$formula, "\n",
    "Data dimensions: ", dim(fmm$data), "\n",
    "Total size of functional domain: ", length(fmm$out_index), "\n",
    "Size of functional domain in fit: ", length(fmm$argvals), "\n"
  )
}

#' Printing the fastFMMconc object
#'
#' Simple method to look at the relevant fields.
#'
#' @param fmm Object to print
#'
#' @return Basic list of entries
#' @method print fastFMM
#' @export

print.fastFMMconc <- function(fmm) {
  paste0(
    "Formula: ", fmm$formula, "\n",
    "Data dimensions: ", dim(fmm$data), "\n",
    "Total size of functional domain: ", length(fmm$out_index), "\n",
    "Size of functional domain in fit: ", length(fmm$argvals), "\n"
  )
}
