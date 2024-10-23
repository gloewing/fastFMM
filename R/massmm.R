#' Generic "massmm" model fit
#'
#' Fits the model over the functional domain based on the object class.
#'
#' @param fmm An object that is or inherits from the "fastFMM" class.
#' @param ... Additional arguments.
#'
#' @return Results depends on the dispatched method.
#' @export

massmm <- function(fmm, ...) {
  UseMethod("massmm")
}

#' Massively univariate non-concurrent fit
#'
#' Fits separate models over the functional domain.
#'
#' @method massmm fastFMM
#' @param fmm Object of "fastFMM" class or its inheritors.
#' @param parallel Boolean, passed from `fui`.
#' @param n_cores Number of cores for parallelization. Defaults to 1.
#'
#' @return A list containing results from the massive univariate step.
#'
#' @importFrom lme4 VarCorr getME
#' @importFrom stats model.matrix
#' @importFrom Matrix t
#' @export

massmm.fastFMM <- function(fmm, parallel, n_cores) {
  # Generate the univariate fits
  mass_list <- massmm_apply(fmm, parallel, n_cores)
  res <- massmm_outs(mass_list, fmm)
  res$analytic <- fmm$analytic

  # If not analytic, stop here
  if (!fmm$analytic)
    return(res)

  var_random <- res$var_random
  se_mat <- suppressMessages(do.call(cbind, lapply(mass_list, '[[', 'se_mat')))
  sigmaesqHat <- var_random["var.Residual", , drop = FALSE]
  # sigmausqHat <- var_random[
  #   which(rownames(var_random) != "var.Residual"), , drop = FALSE
  # ]

  # Return random effects
  randeffs <- NULL
  if (fmm$randeffs) {
    randeffs <- suppressMessages(
      simplify2array(
        lapply(
          lapply(massmm, '[[', 'randeffs'),
          function(x) as.matrix(x[[1]])
        )
      )
    )
  }

  # Fit a sample model at an arbitrary point to obtain a design matrix
  data <- fmm$data
  form <- as.character(fmm$formula)
  out_index <- grep(paste0("^", form[2]), names(data))
  data$Yl <- unclass(data[, out_index][, 1]) # arbitrary index
  fit_uni <- unimm_lmer(
    as.formula(paste0("Yl ~ ", form[3])), data, fmm$family
  )
  designmat <- stats::model.matrix(fit_uni)

  # Names of random effects
  varcorr_df <- as.data.frame(lme4::VarCorr(fit_uni))
  varcorr_df <- varcorr_df[varcorr_df$grp != "Residual", 1:3]
  ztlist <- sapply(lme4::getME(fit_uni, "Ztlist"), Matrix::t)

  # Check if group contains ":" (hierarchical structure) that requires the
  # group to be specified
  group <- mass_list[[1]]$group ## group name in the data
  subj_id <- fmm$subj_id
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
    sigmausqHat = res$sigmausqHat,
    se_mat = se_mat,
    var_random = var_random,
    varcorr_df = varcorr_df,
    randeffs = randeffs,
    ztlist = ztlist,
    designmat = designmat,
    group = group,
    subj_id = subj_id,
    randintercept = randintercept,
    out_index = out_index
  )
  res <- c(res, res_analytic)
  return(res)
}

#' Massively univariate concurrent fit
#'
#' Fits separate models over the functional domain for concurrent models.
#'
#' @method massmm fastFMMconc
#' @param fmm Object of "fastFMM" class or its inheritors.
#' @param parallel Boolean, passed from `fui`.
#' @param n_cores Number of cores for parallelization. Defaults to 1.
#'
#' @return A list of outputs from `unimm`.
#'
#' @importFrom lme4 VarCorr getME
#' @importFrom stats model.matrix
#' @export

massmm.fastFMMconc <- function(fmm, parallel, n_cores) {
  # Generate the univariate fits
  mass_list <- massmm_apply(fmm, parallel, n_cores)
  # print(mass_list)
  res <- massmm_outs(mass_list, fmm)
  res$analytic <- fmm$analytic

  # If not analytic, stop here
  if (!fmm$analytic)
    return(res)

  # If analytic, add outputs for variance estimation
  var_random <- res$var_random
  sigmausqHat <- res$sigmausqHat
  se_mat <- suppressMessages(do.call(cbind, lapply(mass_list, '[[', 'se_mat')))
  sigmaesqHat <- var_random["var.Residual", , drop = FALSE]
  # sigmausqHat <- var_random[
  #   which(rownames(var_random) != "var.Residual"), , drop = FALSE
  # ]

  # Collect random effects
  randeffs <- NULL
  if (fmm$randeffs) {
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
  designmat <- lapply(mass_list, "[[", "designmat")

  # Check if group contains ":" (hierarchical structure) that requires the
  # group to be specified
  group <- mass_list[[1]]$group ## group name in the data
  subj_id <- fmm$subj_id
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
  # AX: There's gotta be an easier way to do this
  # randintercept <- I(
  #   length(fit_uni@cnms) == 1 &
  #     length(fit_uni@cnms[[group]]) == 1 &
  #     fit_uni@cnms[[group]][1] == "(Intercept)"
  # )
  randintercept <- I(
    nrow(mass_list[[1]][["varcorr_df"]]) == 1 &
      mass_list[[1]][["varcorr_df"]][1, "var1"] == "(Intercept)"
  )

  data <- fmm$data
  form <- as.character(fmm$formula)
  out_index <- grep(paste0("^", form[2]), names(data))

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
    out_index = out_index
  )

  res <- c(res, res_analytic)
  return(res)
}

#' Massively univariate parallelization helper
#'
#' A helper function that handles the parallelization of the model fitting.
#' Helps to clean up the code of `massmm` executions. Dispatches `unimm`
#' based on the class of the provided univariate model object.
#'
#' @param fmm Object of "fastFMM" class or its inheritors.
#' @param parallel Boolean, taken from `fui` arguments.
#' @param n_cores Number of cores for parallelization. Defaults to 1.
#'
#' @importFrom parallel mclapply makePSOCKcluster parLapply
#'
#' @return A list of fitted univariate models.
#' @export

massmm_apply <- function(fmm, parallel, n_cores = 1) {
  if (!parallel) {
    res <- lapply(
      fmm$argvals,
      unimm,
      fmm = fmm
    )
    return(res)
  }

  # Windows is not currently compatible with mclapply
  if (.Platform$OS.type != "windows") {
    res <- parallel::mclapply(
      fmm$argvals,
      unimm,
      fmm = fmm,
      mc.cores = n_cores
    )
    return(res)
  }

  # Windows-only parallelization
  # To prevent the cluster from persisting after failure, we use tryCatch
  cl <- parallel::makePSOCKcluster(n_cores)
  res <- tryCatch(
    parLapply(
      cl = cl,
      X = fmm$argvals,
      fun = unimm,
      fmm = fmm
    ),
    warning = function(w) {
      print(paste0("Warning in parallelization of univariate fits:", "\n", w))
    },
    error = function(e) {
      stop(paste0("Error in parallelization of univariate fits:", "\n", e))
      stopCluster(cl)
    },
    finally = {
      # if (!silent) print("Finished fitting univariate models.")
      stopCluster(cl)
    }
  )
  return(res)
}

#' Massively univariate list outputs
#'
#' Helper function that returns basic outputs from `massmm_apply`.
#'
#' @param mass_list The output from `massmm_apply`
#' @param fmm Object of class "fastFMM" with parameters for model fitting.
#'
#' @return A list of betaTilde, the AIC matrix, and residuals (if applicable).
#' @export

massmm_outs <- function(mass_list, fmm) {
  # Obtain betaTilde, fixed effects estimates
  argvals <- fmm$argvals
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
  if (fmm$residuals)
    resids <- suppressMessages(
      do.call(cbind, lapply(mass_list, '[[', 'residuals'))
    )

  # Add smoothed coefficients
  # If analytic, add outputs for variance estimation
  var_random <- t(do.call(rbind, lapply(mass_list, '[[', 'var_random')))
  sigmausqHat <- var_random[
    which(rownames(var_random) != "var.Residual"), , drop = FALSE
  ]

  return(
    list(
      betaTilde = betaTilde,
      AIC_mat = AIC_mat,
      residuals = resids,
      var_random = var_random,
      sigmausqHat = sigmausqHat
    )
  )
}

