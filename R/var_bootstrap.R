#' Analytic variance calculation
#'
#' Helper for `fui`. Bootstrapped variance calculation.
#'
#' @param mum Massively univariate model output of class "massmm"
#' @param umm Univariate mixed model setup of class "unimm"
#' @param nknots_min Integer passed from `fui`.
#' @param nknots_min_cov Integer passed from `fui`.
#' @param betaHat Numeric matrix of smoothed coefficients
#' @param data Data frame of values to fit
#' @param L integer, number of points on functional domain
#' @param num_boots Integer, number of bootstrap replications.
#' @param parallel Logical, whether to use parallel processing
#' @param n_cores Integer, number of cores for parallelization.
#' @param smooth_method Character, passed from `fui`
#' @param splines Character, passed from `fui`
#' @param silent Logical, suppresses messages when `TRUE`. Passed from `fui`.
#'
#' @return List of final outputs of `fui`
#'
#' @import Matrix
#' @importFrom stats smooth.spline quantile
#' @importFrom Rfast spdinv rowMaxs
#' @importFrom parallel mclapply
#' @importFrom methods new
#' @importFrom mvtnorm rmvnorm

var_bootstrap <- function(
  mum,
  umm,
  nknots_min,
  nknots_min_cov,
  betaHat,
  data,
  L,
  num_boots,
  parallel,
  n_cores,
  smooth_method,
  splines,
  silent
) {

  if (!silent) print("Step 3: Inference (Bootstrap)")

  # Check to see if group contains ":" which indicates hierarchical structure
  # and group needs to be specified
  group <- mum$group
  subj_id <- mum$subj_id
  argvals <- mum$argvals

  if (grepl(":", group, fixed = TRUE)) {
    if (is.null(subj_id)) {
      # assumes the ID name is to the right of the ":"
      group <- str_remove(group, ".*:")
    } else if (!is.null(subj_id)) {
      group <- subj_id # use user specified if it exists
    } else {
      message("You must specify the argument: ID")
    }
  }
  ID.number <- unique(data[, group])
  idx_perm <- t(
    replicate(
      num_boots,
      sample.int(length(ID.number), length(ID.number), replace = TRUE)
    )
  )
  B <- num_boots
  betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))

  family <- umm$family
  if (is.null(boot_type)) {
    # default bootstrap type if not specified
    if (family == "gaussian") {
      boot_type <- ifelse(length(ID.number) <= 10, "reb", "wild")
    } else {
      boot_type <- "cluster"
    }
  }

  if (family != "gaussian" & boot_type %in% c("wild", "reb")) {
    stop(
      paste0(
        'Non-gaussian outcomes only supported for some bootstrap procedures.',
        '\n',
        'Set argument `boot_type` to one of the following:', '\n',
        ' "parametric", "semiparametric", "cluster", "case", "residual"'
      )
    )
  }
  message(paste("Bootstrapping Procedure:", as.character(boot_type)))
  # original way
  if (boot_type == "cluster") {
    # Do bootstrap
    pb <- progress_bar$new(total = B)
    for(boots in 1:B) {
      pb$tick()
      # take one of the randomly sampled (and unique) combinations
      sample.ind <- idx_perm[boots,]
      dat.ind <- new_ids <- vector(length = length(sample.ind), "list")
      for(ii in 1:length(sample.ind)) {
        dat.ind[[ii]] <- which(data[, group] == ID.number[sample.ind[ii]])
        # subj_b is now the pseudo_id
        new_ids[[ii]] <- rep(ii, length(dat.ind[[ii]]))
      }
      dat.ind <- do.call(c, dat.ind)
      new_ids <- do.call(c, new_ids)
      df2 <- data[dat.ind, ] # copy dataset with subset of rows we want
      df2[,subj_id] <- new_ids # replace old IDs with new IDs

      fit_boot <- fui(
        formula = umm$formula,
        data = df2,
        family = umm$family,
        argvals = argvals,
        var = FALSE,
        parallel = FALSE,
        silent = TRUE,
        nknots_min = nknots_min,
        nknots_min_cov = nknots_min_cov,
        smooth_method = smooth_method,
        splines = splines,
        residuals = umm$residuals,
        subj_id = mum$subj_id,
        n_cores = n_cores,
        REs = FALSE
      )

      betaHat_boot[, , boots] <- fit_boot$betaHat
    }
    rm(fit_boot, df2, dat.ind, new_ids)
  } else {

    # lmeresampler() way
    # Use original amount. Do not constrain by number of unique resampled
    # types here because we cannot construct rows to resample.
    B <- num_boots
    betaHat_boot <- betaTilde_boot <- array(
      NA, dim = c(nrow(betaHat), ncol(betaHat), B)
    )
    # , as.character(boot_type)
    if (!silent)
      print("Step 3.1: Bootstrap resampling-")

    pb <- progress_bar$new(total = L)
    for(l in 1:L) {
      pb$tick()
      data$Yl <- unclass(data[,out_index][,argvals[l]])
      fit_uni <- suppressMessages(
        lmer(
          formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
          data = data,
          control = lmerControl(
            optimizer = "bobyqa", optCtrl = list(maxfun = 5000))
        )
      )

      # set seed to make sure bootstrap replicate (draws) are correlated
      # across functional domains
      set.seed(seed)

      if (boot_type == "residual") {
        # for residual bootstrap to avoid singularity problems
        boot_sample <- lmeresampler::bootstrap(
          model = fit_uni,
          B = B,
          type = boot_type,
          rbootnoise = 0.0001
        )$replicates
        betaTilde_boot[,l,] <- t(as.matrix(boot_sample[,1:nrow(betaHat)]))
      }else if (boot_type %in% c("wild", "reb", "case") ) {
        # for case
        flist <- lme4::getME(fit_uni, "flist")
        re_names <- names(flist)
        clusters_vec <- c(rev(re_names), ".id")
        # for case bootstrap, we only resample at first (subject level)
        # because doesn't make sense to resample within-cluster for
        # longitudinal data
        resample_vec <- c(TRUE, rep(FALSE, length(clusters_vec) - 1))
        boot_sample <- lmeresampler::bootstrap(
          model = fit_uni,
          B = B,
          type = boot_type,
          resample = resample_vec, # only matters for type = "case"
          hccme = "hc2", # wild bootstrap
          aux.dist = "mammen", # wild bootstrap
          reb_type = 0
        )$replicates # for reb bootstrap only
        betaTilde_boot[, l, ] <- t(as.matrix(boot_sample[,1:nrow(betaHat)]))
      } else {
        use.u <- ifelse(boot_type == "semiparametric", TRUE, FALSE)
        betaTilde_boot[,l,] <- t(
          lme4::bootMer(
            x = fit_uni, FUN = function(.) {fixef(.)},
            nsim = B,
            seed = seed,
            type = boot_type,
            use.u = use.u
          )$t
        )
      }
    }

    suppressWarnings(rm(boot_sample, fit_uni))
    # smooth across functional domain
    if (!silent)   print("Step 3.2: Smooth Bootstrap estimates")
    for(b in 1:B) {
      betaHat_boot[,,b] <- t(
        apply(
          betaTilde_boot[, , b],
          1,
          function(x)
            gam(
              x ~ s(argvals, bs = splines, k = (nknots + 1)),
              method = smooth_method
            )$fitted.values
        )
      )
    }

    rm(betaTilde_boot)
  }

  # Obtain bootstrap variance
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  ## account for within-subject correlation
  for(r in 1:nrow(betaHat)) {
    betaHat.var[,,r] <- 1.2 * var(t(betaHat_boot[r, , ]))
  }

  # Obtain qn to construct joint CI using the fast approach
  if (!silent)
    print("Step 3.3: Joint confidence interval construction")
  qn <- rep(0, length = nrow(betaHat))
  ## sample size in simulation-based approach
  N <- 10000
  # set seed to make sure bootstrap replicate (draws) are correlated across
  # functional domains
  set.seed(seed)
  for(i in 1:length(qn)) {
    est_bs <- t(betaHat_boot[i,,])
    # suppress sqrt(Eigen$values) NaNs
    fit_fpca <- suppressWarnings(
      refund::fpca.face(est_bs, knots = nknots_fpca)
    )
    ## extract estimated eigenfunctions/eigenvalues
    phi <- fit_fpca$efunctions
    lambda <- fit_fpca$evalues
    K <- length(fit_fpca$evalues)

    ## simulate random coefficients
    # generate independent standard normals
    theta <- matrix(stats::rnorm(N*K), nrow=N, ncol=K)
    if (K == 1) {
      # scale to have appropriate variance
      theta <- theta * sqrt(lambda)
      # simulate new functions
      X_new <- tcrossprod(theta, phi)
    } else {
      # scale to have appropriate variance
      theta <- theta %*% diag(sqrt(lambda))
      # simulate new functions
      X_new <- tcrossprod(theta, phi)
    }
    # add back in the mean function
    x_sample <- X_new + t(fit_fpca$mu %o% rep(1,N))
    # standard deviation: apply(x_sample, 2, sd)
    Sigma_sd <- Rfast::colVars(x_sample, std = TRUE, na.rm = FALSE)
    x_mean <- colMeans(est_bs)
    un <- rep(NA, N)
    for(j in 1:N) {
      un[j] <- max(abs((x_sample[j, ] - x_mean) / Sigma_sd))
    }
    qn[i] <- stats::quantile(un, 0.95)
  }

  if (!silent)
    message(
      paste0(
        "Complete!", "\n",
        " - Use plot_fui() function to plot estimates", "\n",
        " - For more information, run the command:  ?plot_fui"
      )
    )

  return(
    list(
      betaHat = betaHat,
      betaHat.var = betaHat.var,
      qn = qn,
      aic = AIC_mat,
      residuals = resids,
      bootstrap_samps = B,
      argvals = argvals
    )
  )
}
