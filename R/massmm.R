#' Generic "massmm" model fit
#'
#' Fits the model over the functional domain. The function `fit_unimm` is
#' dispatched based on the `uni` object class.
#'
#' @param uni_model An object that is or inherits from the "unimm" class.
#' @param ... Additional arguments.
#'
#' @return Results depends on the dispatched method.

massmm <- function(uni_model, ...) {
  UseMethod("massmm")
}

#' Massively univariate non-concurrent fit
#'
#' Fits separate models over the functional domain.
#'
#' @param uni_model Object of "unimm" class or its inheritors.
#' @param argvals Integer vector of locations to fit univariate models.
#' @param data Data frame to fit.
#' @param parallel Boolean, taken from `fui` arguments.
#' @param num_cores Number of cores for parallelization. Defaults to 1.
#'
#' @return A list of outputs from `fit_unimm`.
#'
#' @importFrom lme4 VarCorr getME
#' @importFrom stats model.matrix

massmm.unimm <- function(uni_model, argvals, data, parallel) {
  # Generate the univariate fits
  mass_list <- massmm_apply(uni_model, argvals, data, parallel, num_cores)
  res <- massmm_outs(mass_list, uni_model)

  # If not analytic, stop here
  if (!uni_model$analytic) {
    class(res) <- "massmm"
    return(res)
  }

  # If analytic, add additonal outputs for variance estimation
  var_random <- t(do.call(rbind, lapply(mass_list, '[[', 'var_random')))
  se_mat <- suppressMessages(do.call(cbind, lapply(mass_list, '[[', 'se_mat')))
  sigmaesqHat <- var_random["var.Residual", , drop = FALSE]
  sigmausqHat <- var_random[
    which(rownames(var_random) != "var.Residual"), , drop = FALSE
  ]

  # Collect random effects
  randeffs <- NULL
  if (uni_model$REs) {
    randeffs <- suppressMessages(
      simplify2array(
        lapply(
          lapply(massmm, '[[', 'randeffs'),
          function(x) as.matrix(x[[1]])
        )
      )
    )
  }

  # Fit a fake model to obtain a design matrix
  data$Yl <- unclass(data[, out_index][, 1]) # arbitrary index
  fit_uni <- unimm_lmer(
    as.formula(paste0("Yl ~ ", model_formula[3])), data, uni_model$family
  )
  designmat <- stats::model.matrix(fit_uni)

  # Names of random effects
  varcorr_df <- as.data.frame(lme4::VarCorr(fit_uni))
  varcorr_df <- varcorr_df[varcorr_df$grp != "Residual", 1:3]
  ztlist <- sapply(lme4::getME(fit_uni, "Ztlist"), t)

  # Check if group contains ":" (hierarchical structure) that requires the
  # group to be specified
  group <- mass_list[[1]]$group ## group name in the data
  if (grepl(":", group, fixed = TRUE)) {
    if (is.null(subj_id)) {
      # Assumes the ID name is to the right of the ":"
      group <- sub(".*:", "", group)
    } else {
      # Use user-specified ID if it exists
      group <- subj_id
    }
  }

  if (is.null(subj_id))
    subj_id <- group

  # Condition for using G_estimate_randint
  randintercept <- I(
    length(fit_uni@cnms) == 1 &
      length(fit_uni@cnms[[group]]) == 1 &
      fit_uni@cnms[[group]][1] == "(Intercept)"
  )

  # Check for what needs to be returned
  # This will get passed to smoothing
  res_analytic <- list(
    sigmaesqHat = sigmaesqHat,
    sigmausqHat = sigmausqHat,
    se_mat = se_mat,
    var_random = var_random,
    varcorr_df = varcorr_df,
    randeffs = randeffs,
    ztlist = ztlist,
    designmat = designmat,
    group = group,
    subj_id = subj_id,
    randintercept = randintercept,
    argvals = argvals,
    out_index = out_index
  )
  res <- c(res, res_analytic)
  class(res) <- "massmm"
  return(res)
}

#' Massively univariate concurrent fit
#'
#' Fits separate models over the functional domain for concurrent models.
#'
#' @param uni_model Object of "unimm" class or its inheritors.
#' @param argvals Integer vector of locations to fit univariate models.
#' @param data Data frame to fit.
#' @param parallel Boolean, taken from `fui` arguments.
#' @param num_cores Number of cores for parallelization. Defaults to 1.
#'
#' @return A list of outputs from `fit_unimm`.
#'
#' @importFrom lme4 VarCorr getME
#' @importFrom stats model.matrix

massmm.unimm_conc <- function(uni_model, argvals, data, parallel) {
  # Generate the univariate fits
  mass_list <- massmm_apply(uni_model, argvals, data, parallel, num_cores)
  res <- massmm_outs(mass_list, uni_model)

  # If not analytic, stop here
  if (!uni_model$analytic) {
    class(res) <- "massmm_conc"
    return(res)
  }

  # If analytic, add outputs for variance estimation
  var_random <- t(do.call(rbind, lapply(mass_list, '[[', 'var_random')))
  se_mat <- suppressMessages(do.call(cbind, lapply(mass_list, '[[', 'se_mat')))
  sigmaesqHat <- var_random["var.Residual", , drop = FALSE]
  sigmausqHat <- var_random[
    which(rownames(var_random) != "var.Residual"), , drop = FALSE
  ]

  # Collect random effects
  randeffs <- NULL
  if (uni_model$REs) {
    randeffs <- suppressMessages(
      simplify2array(
        lapply(
          lapply(massmm, '[[', 'randeffs'),
          function(x) as.matrix(x[[1]])
        )
      )
    )
  }

  # Names of random effects
  varcorr_df <- lapply(mass_list, "[[", "varcorr_df")
  ztlist <- lapply(mass_list, "[[", "ztlist")
  ztlist <- lapply(mass_list, "[[", "designmat")

  # Check if group contains ":" (hierarchical structure) that requires the
  # group to be specified
  group <- mass_list[[1]]$group ## group name in the data
  if (grepl(":", group, fixed = TRUE)) {
    if (is.null(subj_id)) {
      # Assumes the ID name is to the right of the ":"
      group <- sub(".*:", "", group)
    } else {
      # Use user-specified ID if it exists
      group <- subj_id
    }
  }

  if (is.null(subj_id))
    subj_id <- group

  # Condition for using G_estimate_randint
  randintercept <- I(
    length(fit_uni@cnms) == 1 &
      length(fit_uni@cnms[[group]]) == 1 &
      fit_uni@cnms[[group]][1] == "(Intercept)"
  )

  # Check for what needs to be returned
  # This will get passed to smoothing
  res_analytic <- list(
    sigmaesqHat = sigmaesqHat,
    sigmausqHat = sigmausqHat,
    se_mat = se_mat,
    var_random = var_random,
    varcorr_df = varcorr_df,
    randeffs = randeffs,
    ztlist = ztlist,
    designmat = designmat,
    group = group,
    subj_id = subj_id,
    randintercept = randintercept,
    argvals = argvals,
    out_index = out_index
  )
  # Manually set this to false due to lack of shortcut coded
  res_analytic$randintercept <- FALSE
  res <- c(res, res_analytic)
  class(res) <- "massmm_conc"
  return(res)
}

#' Massively univariate parallelization helper
#'
#' A helper function that handles the parallelization of the model fitting.
#' Helps to clean up the code of `fit_massmm` executions. Dispatches `fit_unimm`
#' based on the class of the provided univariate model object.
#'
#' @param uni_model Object of "unimm" class or its inheritors.
#' @param argvals Numeric vector of locations to fit univariate models.
#' @param data Data frame to fit.
#' @param parallel Boolean, taken from `fui` arguments.
#' @param num_cores Number of cores for parallelization. Defaults to 1.
#'
#' @importFrom parallel mclapply makePSOCKcluster parLapply
#'
#' @return A list of fitted univariate models.

massmm_apply <- function(uni_model, argvals, data, parallel, num_cores = 1) {
  if (!parallel) {
    res <- lapply(
      argvals,
      fit_unimm,
      data = data,
      obj = uni_model
    )
    return(res)
  }

  # Windows is not currently compatible with mclapply
  if (.Platform$OS.type != "windows") {
    res <- parallel::mclapply(
      argvals,
      fit_unimm,
      data = data,
      obj = uni_model,
      mc.cores = num_cores
    )
    return(res)
  }

  # Windows-only parallelization
  # To prevent the cluster from persisting after failure, we use tryCatch
  cl <- parallel::makePSOCKcluster(num_cores)
  res <- tryCatch(
    parLapply(
      cl = cl,
      X = argvals,
      fun = fit_unimm,
      data = data,
      obj = uni_model
    ),
    warning = function(w) {
      print(paste0("Warning in parallelization of univariate fits:", "\n", w))
    },
    error = function(e) {
      stop(paste0("Error in parallelization of univariate fits:", "\n", e))
      stopCluster(cl)
    },
    finally = {
      if (!silent) print("Finished fitting univariate models.")
      stopCluster(cl)
    }
  )
  return(res)
}

#' Massively univariate list outputs
#'
#' Helper function that returns basic outputs from `massmm_apply`.
#'
#' @param massmm_list The output from `massmm_apply`
#' @param uni_model A "unimm" object with model parameters.
#'
#' @return A list of betaTilde, the AIC matrix, and residuals (if applicable).

massmm_outs <- function(massmm_list, uni_model) {
  # Obtain betaTilde, fixed effects estimates
  betaTilde <- t(do.call(rbind, lapply(mass_list, '[[', 1)))
  colnames(betaTilde) <- argvals

  # Obtain AIC and BIC
  mod_aic <- do.call(c, lapply(mass_list, '[[', 'aic'))
  mod_bic <- do.call(c, lapply(mass_list, '[[', 'bic'))
  mod_caic <- do.call(c, lapply(mass_list, '[[', 'caic'))
  AIC_mat <- cbind(mod_aic, mod_bic, mod_caic)
  colnames(AIC_mat) <- c("AIC", "BIC", "cAIC")

  # Store residuals if resids == TRUE
  resids <- NA
  if (uni_model$residuals)
    resids <- suppressMessages(
      do.call(cbind, lapply(mass_list, '[[', 'residuals'))
    )

  return(
    list(
      betaTilde = betaTilde,
      AIC_mat = AIC_mat,
      resids = resids,
      analytic = uni_model$analytic
    )
  )
}

