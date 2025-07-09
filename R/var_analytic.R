#' Analytic variance calculation
#'
#' Helper for `fui`. Analytic calculation for Gaussian models.
#'
#' @param fmm Object of class "fastFMM".
#' @param mum Massively univariate model output of class "massmm"
#' @param smoothed Outputs from the smoothing Step 2
#' @param MoM integer, indicates type of MoM estimator (1 or 2)
#' @param non_neg Integer, indicates type of non-negativity calculation
#' @param nknots_cov Integer, number of knots for splines
#' @param seed Integer random seed
#' @param parallel Logical, whether to use parallel processing
#' @param n_cores Integer number of cores for parallelization
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
#' @export

var_analytic <- function(
  fmm,
  mum,
  smoothed,
  MoM,
  non_neg,
  nknots_cov,
  seed,
  parallel,
  n_cores,
  silent
) {

  # 0 Preparation of H, R #####################################################

  # Variance estimates of random components: H(s), R(s)
  if (!silent) print("Step 3: Inference (Analytic)")
  if (!silent) print("Step 3.0: Preparation")

  # Source variables
  betaHat <- smoothed$betaHat
  HHat <- smoothed$HHat
  lambda <- smoothed$lambda
  S <- smoothed$S
  B <- smoothed$B

  # AX: add imputation here
  data <- fmm$data
  randintercept <- mum$randintercept
  argvals <- fmm$argvals
  designmat <- mum$designmat
  p <- nrow(mum$betaTilde)
  L <- length(fmm$argvals)
  RHat <- t(
    apply(
      mum$sigmaesqHat, 1,
      function(b) stats::smooth.spline(x = argvals, y = b)$y
    )
  )
  RHat[which(RHat < 0)] <- 0

  # 1 MoM estimator of G #######################################################

  if (randintercept) {

    # 1.1 Preparation A: Random intercept case =================================
    # Includes potential NNLS correction for diagonals (for variance terms)

    GTilde <- G_estimate_randint(fmm, designmat, betaHat, silent)

    if (!silent)
      print("Step 3.1.1: Preparation A")

    diag(GTilde) <- HHat[1, ] # L x L matrix
    # Fast bivariate smoother
    # nknots_min

    if (!silent) print("Step 3.1.2: Smooth G")

    GHat <- fbps_cov(
      GTilde,
      search.length = 100,
      knots = nknots_cov
    )$Yhat
    diag(GHat)[which(diag(GHat) < 0)] <- diag(GTilde)[which(diag(GHat) < 0)]

  } else {

    # 1.1 Preparation B ==========================================================

    ## Runs if there are more random effects than the intercept
    if (!silent)
      print("Step 3.1.1: Preparation B")

    # generate design matrix for G(s_1, s_2)
    # Method of Moments Linear Regression calculation

    ## Method of Moments estimator for G() with potential NNLS correction
    # for diagonals (for variance terms)

    G_est_out <- G_estimate(fmm, mum, betaHat, HHat, non_neg, MoM, silent)
    GTilde <- G_est_out$GTilde
    data_cov <- G_est_out$data_cov

    # 1.2 Smooth G =============================================================

    if (!silent) print("Step 3.1.2: Smooth G")

    ## smooth GHat
    GHat <- GTilde
    for (r in 1:nrow(HHat)) {
      diag(GTilde[r, , ]) <- HHat[r, ]
      GHat[r, , ] <- fbps_cov(
        GTilde[r, , ],
        search.length = 100,
        knots = nknots_cov # nknots_min
      )$Yhat
    }
    diag(GHat[r, , ])[which(diag(GHat[r, , ]) < 0)] <-
      diag(GTilde[r, , ])[which(diag(GHat[r, , ]) < 0)]
  }

  # 2 First step ###############################################################

  if (!silent) print("Step 3.2: First step")

  # Calculate the intra-location variance of betaTilde: Var(betaTilde(s))

  # Obtain the corresponding rows of each subject
  obs_ind <- list()
  group <- mum$group
  id_list <- unique(data[, group])
  for(id in id_list) {
    obs_ind[[as.character(id)]] <- which(data[, group] == id)
  }

  # Concatenate vector of 1s to Z because used that way below
  HHat_trim <- NA

  if (!randintercept) {
    Z <- data_cov$Z_orig
    if (is.null(Z)) Z <- data_cov$Z # Concurrent case
    qq <- ncol(Z)
    HHat_trim <- array(NA, c(qq, qq, L)) # array for Hhat
  }

  # Create betaTilde_theo_var to store variance of betaTilde
  betaTilde_theo_var <- array(NA, dim = c(p, p, L))
  # Create XTVinvZ_all to store all XTVinvZ used in the covariance calculation
  XTVinvZ_all <- vector(length = L, "list")
  ## arbitrarily start find indices
  resStart <- cov_organize_start(HHat[,1])
  res_template <- resStart$v_list_template # index template
  template_cols <- ncol(res_template)

  # Parallelize the calculation
  if (parallel) {
    # Windows is incompatible with mclapply
    if(.Platform$OS.type == "windows") {
      cl <- parallel::makePSOCKcluster(n_cores)
      # Provide necessary objects to the cluster
      clusterExport(cl, "eigenval_trim", environment())

      # tryCatch ensures the cluster is stopped even during errors
      massVar <- tryCatch(
        parLapply(
          cl = cl,
          X = argvals,
          fun = var_parallel,
          fmm = fmm,
          mum = mum,
          Z = Z,
          RHat = RHat,
          HHat = HHat,
          id_list = id_list,
          obs_ind = obs_ind,
          res_template = res_template
        ),
        warning = function(w) {
          print(
            paste0(
              "Warning when parallelizing variance calculation:", "\n", w
            )
          )
        },
        error = function(e) {
          stop(
            paste0("Error when parallelizing variance calculations:", "\n", e)
          )
          stopCluster(cl)
        },
        finally = {
          stopCluster(cl)
          if (!silent)
            print("Finished parallelized variance calculation.")

        }
      )
    } else {
      # if not Windows, use mclapply
      massVar <- parallel::mclapply(
        X = argvals,
        FUN = var_parallel,
        mc.cores = n_cores,
        fmm = fmm,
        mum = mum,
        Z = Z,
        RHat = RHat,
        HHat = HHat,
        id_list = id_list,
        obs_ind = obs_ind,
        res_template = res_template
      )
    }
  } else {
    massVar <- lapply(
      argvals,
      var_parallel,
      fmm = fmm,
      mum = mum,
      Z = Z,
      RHat = RHat,
      HHat = HHat,
      id_list = id_list,
      obs_ind = obs_ind,
      res_template = res_template
    )
  }

  XTVinvZ_all <- lapply(argvals, function(s)
    massVar[[s]]$XTVinv[lengths(massVar[[s]]$XTVinv) != 0])
  betaTilde_theo_var <- lapply(argvals, function(s)
    massVar[[s]]$betaTilde_theo_var)
  betaTilde_theo_var <- array(
    simplify2array(betaTilde_theo_var), dim = c(p, p, L)
  )

  suppressWarnings(rm(massVar, resStart, res_template, template_cols))

  # Calculate the inter-location covariance of betaTilde:
  # Cov(betaTilde(s_1), betaTilde(s_2))
  ## Create betaTilde_theo_cov to store covariance of betaTilde
  betaTilde_theo_cov <- array(NA, dim = c(p, p, L, L))
  if (randintercept) {
    resStart <- cov_organize_start(GHat[1, 2]) # arbitrarily start
  } else {
    resStart <- cov_organize_start(GHat[, 1, 2]) # arbitrarily start
  }

  # 2.1 First step =============================================================

  if (!silent) print("Step 3.2.1: First step")

  res_template <- resStart$v_list_template # index template
  template_cols <- ncol(res_template)
  ## Calculate Cov(betaTilde) for each pair of location
  for (i in 1:L) {
    for (j in i:L) {
      tmp <- matrix(0, nrow = p, ncol = p) ## store intermediate part
      if (randintercept) {
        G_use <- GHat[i, j]
      } else {
        G_use <- eigenval_trim(
          matrix(c(GHat[, i, j], 0)[res_template], template_cols)
        )
      }

      for (id in id_list) {
        tmp <- tmp +
          XTVinvZ_all[[i]][[as.character(id)]] %*%
          tcrossprod(G_use, XTVinvZ_all[[j]][[as.character(id)]])
      }

      ## Calculate covariance using XTVinvX and tmp to save memory
      betaTilde_theo_cov[, , i, j] <- betaTilde_theo_var[, , i] %*%
        tmp %*%
        betaTilde_theo_var[, , j]
    }
  }
  suppressWarnings(
    rm(
      Z,
      XTVinvZ_all,
      resStart,
      res_template,
      template_cols
    )
  )

  # 3 Second step ##############################################################

  if (!silent) print("Step 3.3: Second step")

  # Intermediate step for covariance estimate
  betaTilde_var_s <- array(NA, dim = c(L, L, p))
  for(j in 1:p) {
    for(r in 1:L) {
      for(t in 1:L) {
        if (t == r) {
          betaTilde_var_s[r, t, j] <- betaTilde_theo_var[j, j, r]
        } else {
          betaTilde_var_s[r, t, j] <- betaTilde_theo_cov[
            j, j, min(r, t), max(r, t)
          ]
        }
      }
    }
  }

  # Calculate the inter-location covariance of betaHat:
  # Cov(betaHat(s_1), betaHat(s_2))
  betaHat_var <- array(NA, dim = c(L, L, p))
  for(r in 1:p) {
    M <- B %*%
      Matrix::tcrossprod(
        solve(Matrix::crossprod(B) + lambda[r] * S), B
      ) +
      matrix(1 / L, nrow = L, ncol = L)
    raw_var <- M %*% Matrix::tcrossprod(betaTilde_var_s[, , r], M)
    ## trim eigenvalues to make final variance matrix PSD
    betaHat_var[, , r] <- eigenval_trim(raw_var)
  }

  # Obtain qn to construct joint CI
  qn <- rep(0, length = nrow(betaHat))
  N <- 10000 ## sample size in simulation-based approach
  zero_vec <- rep(0, length(betaHat[1,]))
  set.seed(seed)

  for(i in 1:length(qn)) {
    Sigma <- betaHat_var[, , i]
    sqrt_Sigma <- sqrt(diag(Sigma))
    S_scl <- Matrix::Diagonal(x = 1 / sqrt_Sigma)
    Sigma <- as.matrix(S_scl %*% Sigma %*% S_scl)
    # x_sample <- abs(FastGP::rcpp_rmvnorm_stable(N, Sigma, zero_vec))
    x_sample <- abs(mvtnorm::rmvnorm(N, zero_vec, Sigma))
    un <- Rfast::rowMaxs(x_sample, value = TRUE)
    qn[i] <- stats::quantile(un, 0.95)
  }

  # Decide whether to return design matrix or just set it to NULL
  # if (!design_mat) designmat <- NA
  if (!silent)
    message(
      paste0(
        "Complete!", "\n",
        " - Use plot_fui() function to plot estimates.", "\n",
        " - For more information, run the command:  ?plot_fui"
      )
    )

  return(
    list(
      betaHat = betaHat,
      betaHat_var = betaHat_var,
      qn = qn,
      aic = mum$AIC_mat,
      betaTilde = mum$betaTilde,
      var_random = mum$var_random,
      designmat = mum$designmat,
      residuals = mum$resids,
      H = HHat_trim,
      R = RHat,
      G = GTilde,
      GHat = GHat,
      Z = mum$ztlist,
      argvals = fmm$argvals,
      randeffs = mum$randeffs,
      se_mat = mum$se_mat
    )
  )
}

# Helper functions #############################################################

# Functions related to G generation/estimation are in separate files for
# convenience.

#' select_knots.R from refund package
#'
#' Copied from [select_knots](https://rdrr.io/cran/refund/src/R/select_knots.R)
#' because the original is not exported for use.
#'
#' @param t Numeric
#' @param knots Numeric scalar or vector, the number/numbers of  knots or the
#' vector/vectors of knots for each dimension. Default = 10
#' @param p Numeric, the degrees of B-splines. Default = 3.
#' @param option Character, knot spacing, can be `"equally-spaced"` or
#' `"quantile"`
#'
#' @return Vector of specified knots
#'
#' @importFrom stats quantile
#' @export

select_knots <- function(
    t, knots = 10, p = 3, option = "equally-spaced"
) {
  qs <- seq(0, 1, length = knots + 1)

  if (option == "equally-spaced")
    knots <- (max(t) - min(t)) * qs + min(t)

  if (option == "quantile")
    knots <- as.vector(stats::quantile(t,qs))


  K <- length(knots)
  knots_left <- 2 * knots[1] - knots[p:1 + 1]
  knots_right <- 2 * knots[K] - knots[K - (1:p)]

  return(c(knots_left, knots, knots_right))
}

#' Fast block diagonal generator, taken from Matrix package examples
#'
#' Copyright (C) 2016 Martin Maechler, ETH Zurich. Copied here for more
#' convenient integration.
#'
#' @param lmat Matrix
#'
#' @return Block diagonal version of input
#'
#' @importFrom methods new
#' @export

bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if (!length(lmat))
    return(methods::new("dgCMatrix"))

  stopifnot(
    is.list(lmat),
    is.matrix(lmat[[1]]),
    (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
    all(vapply(lmat, dim, integer(2)) == k)
  ) # all of them

  N <- length(lmat)
  if (N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M = ", N * k)
  M <- as.integer(N * k)

  ## result: an   M x M  matrix of class dgCMatrix
  methods::new(
    "dgCMatrix",
    Dim = c(M,M),
    # 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
    i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
    p = k * 0L:M,
    x = as.double(unlist(lmat, recursive = FALSE, use.names = FALSE))
  )
}

#' Organize covariance matrices
#'
#' Read the data frame of variance-covariance interactiosn and
#' organize the covariance matrices correctly.
#' Finds indices to feed into `cov_organize` function above
#'
#' @param cov_vec Character vector
#'
#' @return List of relevant outputs
#' @export

cov_organize_start <- function(cov_vec) {
  # Assume each set of cov for 2 preceeding variance terms of random effects
  # AX: Check for applicability to more complex random effect interactions.
  # the order returned is the order given by cov_vec (we use (sort below))

  # names of the terms to differentiate variances and  covariances
  nm <- names(cov_vec)
  # find terms that are covariances
  cov_idx <- which( grepl("cov", nm, fixed = TRUE) )
  covs <- length(cov_idx) # number of covariances

  # organize based on number of covariance terms
  if (covs == 0) {
    # if no covariance terms (just simple diagonal matrix)
    var_nm <- groupings <- groups_t <- g_idx_list <- v_list <- NULL
    v_list_template <- diag(1:length(cov_vec))
    # replace 0s with corresponding index
    v_list_template <- apply(
      v_list_template,
      1,
      function(x)
        ifelse(x == 0, max(v_list_template) + 1, x)
    )

  } else {
    # mix of covariance terms and non-covariance terms
    # variance terms
    var_nm <- nm[
      sapply(
        nm,
        function(x)
          unlist(strsplit(x, split='.', fixed=TRUE))[1] == "var")
    ]
    groupings <- unique(
      sapply(
        var_nm,
        function(x)
          unlist(strsplit(x, split='.', fixed=TRUE))[2]
      )
    )
    g_idx_list <- mat_lst <- vector(length = length(groupings), "list")
    cnt <- 0
    # grouping variable for each name
    groups_t <- sapply(
      var_nm,
      function(x)
        unlist(strsplit(x, split='.', fixed=TRUE))[2]
    )
    v_list <- vector(length = length(groupings), "list")

    # iterate through groupings and make sub-matrices
    for(trm in groupings) {
      cnt <- cnt + 1
      # find current grouping (e.g., "id" or "id:session")
      # current grouping
      g_idx <- g_idx_list[[cnt]] <- names(which( groups_t == trm  ))
      # if this is not just a variance term (i.e., cov terms too)
      if (length(g_idx) > 1) {
        # make diagonal matrix with variance terms
        m <- diag(cov_vec[g_idx])
        v_list[[cnt]] <- matrix(NA, ncol = ncol(m), nrow = ncol(m))

        # name after individual random effect variables (covariates)
        trm_var_name <- sapply(
          g_idx,
          function(x)
            unlist(strsplit(x, split='.', fixed=TRUE))[3]
        )
        nm1 <- paste0("cov.", trm, ".")
        for(ll in 1:ncol(m)) {
          for(jj in seq(1,ncol(m) )[-ll] ) {
            # names that are covariance between relevant variables
            v1 <- which(
              nm == paste0(nm1, trm_var_name[ll], ".", trm_var_name[jj])
            )
            v2 <- which(
              nm == paste0(nm1, trm_var_name[jj], ".", trm_var_name[ll])
            )
            v_list[[cnt]][ll,jj] <- c(v1, v2)
          }
        }
        v_list[[cnt]][is.na(v_list[[cnt]])] <- sapply(
          g_idx,
          function(x) which(nm == x)
        )
      }
    }

    # diagonal matrix of indices for template
    v_list_template <- bdiag_m(v_list)
    # replace 0s with corresponding index (just highest + 1)
    v_list_template <- apply(
      v_list_template,
      1,
      function(x) ifelse(x == 0, max(v_list_template) + 1, x)
    )
  }

  return(
    list(
      nm = nm,
      cov_idx = cov_idx,
      covs = covs,
      groupings = groupings,
      g_idx_list = g_idx_list,
      v_list = v_list,
      v_list_template = v_list_template
    )
  )
}

#' Trim eigenvalues and return a positive semidefinite matrix
#'
#' Removes the negative eigenvalues of the given matrix V, returning a positive
#' semidefinite (PSD) matrix.
#'
#' @param V Matrix of dimension n x n
#'
#' @return A Matrix with negative eigenvalues removed
#' @export

eigenval_trim <- function(V) {
  ## trim non-positive eigenvalues to ensure positive semidefinite
  edcomp <- eigen(V, symmetric = TRUE)
  eigen.positive <- which(edcomp$values > 0)
  q <- ncol(V)

  if (length(eigen.positive) == q) {
    # nothing needed here because matrix is already PSD
    return(V)
  } else if (length(eigen.positive) == 0) {
    return(Matrix::tcrossprod(edcomp$vectors[, 1]) * edcomp$values[1])
  } else if (length(eigen.positive) == 1) {
    return(
      Matrix::tcrossprod(as.vector(edcomp$vectors[,1])) *
        as.numeric(edcomp$values[1])
    )
  } else {
    # sum of outer products of eigenvectors, scaled by eigenvalues for all
    # positive eigenvalues
    return(
      matrix(
        edcomp$vectors[, eigen.positive] %*%
          Matrix::tcrossprod(
            diag(edcomp$values[eigen.positive]),
            edcomp$vectors[,eigen.positive]
          ),
        ncol = q
      )
    )
  }
}

#' Copy of pspline.setting.R from refund
#'
#' A slightly modified copy of [pspline.setting](https://rdrr.io/cran/refund/src/R/pspline.setting.R)
#' from `refund`. Copied here because the original function is not exported from
#' the package. Called within `fbps_cov` within `var_analytic`.
#'
#' @param x Marginal data points
#' @param knots The list of interior knots of the numbers of interior knots
#' @param p Degrees for B-splines, default = 3
#' @param m Orders of difference penalty, default = 2
#'
#' @importFrom splines spline.des
#' @import Matrix
#' @export

pspline_setting <- function(
    x,
    knots = select_knots(x, 35),
    p = 3,
    m = 2,
    periodicity = FALSE,
    weight = NULL
) {
  ### design matrix
  K <- length(knots) - 2 * p - 1
  B <- splines::spline.des(
    knots = knots,
    x = x,
    ord = p + 1,
    outer.ok = TRUE
  )$design

  if (periodicity) {
    Bint <- B[, -c(1:p, K + 1:p)]
    Bleft <- B[, 1:p]
    Bright <- B[, K + 1:p]
    B <- cbind(Bint, Bleft + Bright)
  }

  # parameter  m: difference order
  # parameter  p: degree of B-splines
  # parameter  K: number of interior knots
  penalize_difference <-function(m, p, K, periodicity = periodicity){
    c <- rep(0, m + 1)

    for(i in 0:m)
      c[i + 1] = (-1)^(i + 1) * factorial(m) /
        (factorial(i) * factorial(m - i))

    if (!periodicity){
      M <- matrix(0, nrow = K + p - m, ncol = K + p)
      for(i in 1:(K + p - m)) M[i, i:(i + m)] <- c
    }

    if (periodicity){
      M <- matrix(0,nrow = K, ncol = K)
      for(i in 1:(K - m)) M[i, i:(i + m)] <- c
      for(i in (K - m + 1):K) M[i, c(i:K, 1:(m - K + i))] <- c
    }

    return(M)
  }

  P <- penalize_difference(m, p, K, periodicity)
  P1 <- Matrix(P)
  P2 <- Matrix(t(P))
  P <- P2 %*% P1

  MM <- function(A, s, option = 1){
    if (option == 2)
      return(A * (s %*% t(rep(1, dim(A)[2]))))
    if(option == 1)
      return(A * (rep(1, dim(A)[1]) %*% t(s)))
  }

  if(is.null(weight)) weight <- rep(1,length(x))

  B1 <- Matrix(MM(t(B), weight))
  B <- Matrix(B)
  Sig <- B1 %*% B
  eSig <- eigen(Sig)
  V <- eSig$vectors
  E <- eSig$values

  if (min(E) <= 0.0000001)
    E <- E + 0.000001

  Sigi_sqrt = MM(V, 1 / sqrt(E)) %*% t(V)
  #Sigi = V%*%diag(1/E)%*%t(V)

  tUPU <- Sigi_sqrt %*% (P %*% Sigi_sqrt)
  Esig <- eigen(tUPU, symmetric = TRUE)
  U <- Esig$vectors
  s <- Esig$values
  if(!periodicity) s[(K + p - m + 1):(K + p)] <- 0
  if(periodicity) s[K] <- 0
  A <- B %*% (Sigi_sqrt %*% U)

  return(
    list(
      A = A,
      B = B,
      s = s,
      Sigi_sqrt = Sigi_sqrt,
      U = U,
      P = P
    )
  )
}

#' Modified refund::fbps for covariance matrices
#'
#' Altered to speed up `refund::fbps`. Takes advantage of the symmetry of the
#' covariance matrix, i.e., lambda1 == lambda2
#'
#' @param data Data frame of values to fit
#' @param subj Character vector of subject IDs, defaults to `NULL`
#' @param covariates List of data points for each dimension; defaults to `NULL`
#' @param knots Numeric or list of number of knots or the vector of
#' knots for each dimension; defaults to 35.
#' @param knots.option Character knot selection method; defaults to
#' `"equally-spaced"`
#' @param p Numeric degrees of B-splines
#' @param m Numeric order of difference penalty
#' @param lambda User-selected smoothing parameters; defaults to NULL
#' @param selection Character selection of smoothing parameter; defaults to
#' "GCV:
#' @param search.grid boolean; defaults to `TRUE`. If false, uses `optim`.
#' @param search.length Integer no. of equidistant (log) smoothing parameters;
#' defaults to 100.
#' @param method See `stats::optim`; defaults to `"L-BFGS-B"`
#' @param lower,upper Numeric bounds for log smoothing parameter; defaults to
#' `-20, 20`. Passed to `stats::optim`.
#' @param control See `optim`.
#'
#' @return A smoothed matrix
#' @export

fbps_cov <- function(
    data,
    subj = NULL,
    covariates = NULL,
    knots = 35,
    knots.option = "equally-spaced",
    periodicity = c(FALSE,FALSE),
    p = 3,
    m = 2,
    lambda = NULL,
    selection = "GCV",
    search.grid = T,
    search.length = 100,
    method="L-BFGS-B",
    lower = -20,
    upper = 20,
    control = NULL
) {

  ## data dimension
  data_dim <- dim(data)
  n1 <- data_dim[1]
  n2 <- data_dim[2]

  ## subject ID
  if (is.null(subj)) subj <- 1:n2
  subj_unique <- unique(subj)
  I <- length(subj_unique)
  ## covariates for the two axis
  if (!is.list(covariates)) {
    ## if NULL, assume equally distributed
    x <- (1:n1) / n1 - 1 / 2 / n1
    z <- (1:n2) / n2 - 1 / 2 / n2
  }
  if (is.list(covariates)) {

    x <- covariates[[1]]
    z <- covariates[[2]]
  }

  ## B-spline basis setting
  p1 <- rep(p, 2)[1]
  p2 <- rep(p, 2)[2]
  m1 <- rep(m, 2)[1]
  m2 <- rep(m, 2)[2]

  ## knots
  if (!is.list(knots)) {
    K1 <- rep(knots, 2)[1]
    xknots <- select_knots(x, knots = K1, option = knots.option)
    K2 <- rep(knots, 2)[2]
    zknots <- select_knots(z, knots = K2, option = knots.option)
  }

  if (is.list(knots)) {

    xknots <- knots[[1]]
    K1 <- length(xknots) - 1
    knots_left <- 2 * xknots[1] - xknots[p1:1 + 1]
    knots_right <- 2 * xknots[K1] - xknots[K1 - (1:p1)]
    xknots <- c(knots_left, xknots, knots_right)

    zknots<- knots[[2]]
    K2 <- length(zknots)-1
    knots_left <- 2*zknots[1]- zknots[p2:1+1]
    knots_right <- 2*zknots[K2] - zknots[K2-(1:p2)]
    zknots <- c(knots_left,zknots,knots_right)
  }

  Y = data

  # 1 Precalculation for fbps smoothing ########################################

  List <- pspline_setting(x, xknots, p1, m1, periodicity[1])
  A1 <- List$A
  B1 <- List$B
  Bt1 <- Matrix::Matrix(t(as.matrix(B1)))
  s1 <- List$s
  Sigi1_sqrt <- List$Sigi_sqrt
  U1 <- List$U
  A01 <- Sigi1_sqrt%*%U1
  c1 <- length(s1)

  List <- pspline_setting(z, zknots, p2, m2, periodicity[2])
  A2 <- List$A
  B2 <- List$B
  Bt2 <- Matrix(t(as.matrix(B2)))
  s2 <- List$s
  Sigi2_sqrt <- List$Sigi_sqrt
  U2 <- List$U
  A02 <- Sigi2_sqrt%*%U2
  c2 <- length(s2)

  # 2 Select optimal penalty ###################################################

  # Trace of square matrix
  tr <-function(A)
    return(sum(diag(A)))

  Ytilde <- Bt1 %*% (Y %*% B2)
  Ytilde <- t(A01) %*% Ytilde %*% A02
  Y_sum <- sum(Y^2)
  ytilde <- as.vector(Ytilde)

  if (selection=="iGCV") {
    KH <- function(A,B) {
      C <- matrix(0, dim(A)[1], dim(A)[2] * dim(B)[2])
      for(i in 1:dim(A)[1])
        C[i, ] <- Matrix::kronecker(A[i, ], B[i, ])
      return(C)
    }
    G <- rep(0, I)
    Ybar <- matrix(0, c1, I)
    C <- matrix(0, c2, I)

    Y2 <- Bt1%*%Y
    Y2 <- matrix(t(A01) %*% Y2, c1, n2)

    for(i in 1:I) {
      sel <- (1:n2)[subj == subj_unique[i]]
      len <- length(sel)
      G[i] <- len
      Ybar[, i] <- as.vector(matrix(Y2[, sel], ncol = len) %*% rep(1, len))
      C[, i] <- as.vector(t(matrix(A2[sel, ], nrow=len)) %*% rep(1, len))
    }

    g1 <- diag(Ybar %*% diag(G) %*% t(Ybar))
    g2 <- diag(Ybar %*% t(Ybar))
    g3 <- ytilde * as.vector(Ybar %*% diag(G) %*% t(C))
    g4 <- ytilde * as.vector(Ybar %*% t(C))
    g5 <- diag(C %*% diag(G) %*% t(C))
    g6 <- diag(C %*% t(C))
    #cat("Processing completed\n")
  }

  fbps_gcv =function(x) {

    lambda <- exp(x)
    ## two lambda's are the same
    if (length(lambda)==1)
    {
      lambda1 <- lambda
      lambda2 <- lambda
    }
    ## two lambda's are different
    if (length(lambda)==2) {
      lambda1 <- lambda[1]
      lambda2 <- lambda[2]
    }

    sigma2 <- 1 / (1 + lambda2 * s2)
    sigma1 <- 1 / (1 + lambda1 * s1)
    sigma2_sum <- sum(sigma2)
    sigma1_sum <- sum(sigma1)
    sigma <- Matrix::kronecker(sigma2, sigma1)
    sigma.2 <- Matrix::kronecker(sqrt(sigma2), sqrt(sigma1))

    if (selection=="iGCV") {
      d <- 1 / (1 - (1 + lambda1 * s1) / sigma2_sum * n2)
      gcv <- sum((ytilde*sigma)^2) - 2 * sum((ytilde * sigma.2)^2)
      gcv <- gcv + sum(d^2 * g1) - 2 * sum(d * g2)
      gcv <- gcv - 2 * sum(g3 * Matrix::kronecker(sigma2, sigma1 * d^2))
      gcv <- gcv + 4 * sum(g4 * Matrix::kronecker(sigma2, sigma1 * d))
      gcv <- gcv +
        sum(ytilde^2 * Matrix::kronecker(sigma2^2 * g5, sigma1^2 * d^2))
      gcv <- gcv -
        2 * sum(ytilde^2 * Matrix::kronecker(sigma2^2 * g6, sigma1^2 * d))
    }

    if (selection=="GCV") {
      gcv <- sum((ytilde * sigma)^2) - 2 * sum((ytilde * sigma.2)^2)
      gcv <- Y_sum + gcv
      trc <- sigma2_sum*sigma1_sum
      gcv <- gcv / (1 - trc / (n1 * n2))^2
    }
    return(gcv)
  }

  fbps_est <- function(x) {

    lambda <- exp(x)
    ## two lambda's are the same
    if (length(lambda) == 1){
      lambda1 <- lambda
      lambda2 <- lambda
    }
    ## two lambda's are different
    if (length(lambda) ==2){
      lambda1 <- lambda[1]
      lambda2 <- lambda[2]
    }

    sigma2 <- 1 / (1 + lambda2 * s2)
    sigma1 <- 1 / (1 + lambda1 * s1)
    sigma2_sum <- sum(sigma2)
    sigma1_sum <- sum(sigma1)
    sigma <- Matrix::kronecker(sigma2, sigma1)
    sigma.2 <- Matrix::kronecker(sqrt(sigma2), sqrt(sigma1))

    Theta <- A01 %*% diag(sigma1) %*% Ytilde
    Theta <- as.matrix(Theta %*% diag(sigma2) %*% t(A02))
    Yhat <- as.matrix(as.matrix(B1 %*% Theta) %*% Bt2)

    result <- list(
      lambda = c(lambda1, lambda2),
      Yhat = Yhat,
      Theta = Theta,
      setting = list(
        x = list(
          knots = xknots,
          p = p1,
          m = m1
        ),
        z = list(
          knots = zknots,
          p = p2,
          m = m2
        )
      )
    )
    class(result) <- "fbps"
    return(result)
  }

  if (is.null(lambda)) {

    if (search.grid == T) {
      lower2 <- lower1 <- lower[1]
      upper2 <- upper1 <- upper[1]
      search.length2 <- search.length1 <- search.length[1]
      if (length(lower) == 2) lower2 <- lower[2]
      if (length(upper) == 2) upper2 <- upper[2]
      if (length(search.length) == 2) search.length2 <- search.length[2]

      Lambda1 <- seq(lower1,upper1,length = search.length1)
      lambda.length1 <- length(Lambda1)

      GCV = rep(0, lambda.length1)
      for(j in 1:lambda.length1) {
        GCV[j] <- fbps_gcv(c(Lambda1[j], Lambda1[j]))
      }

      location <- which.min(GCV)[1]
      j0 <- location%%lambda.length1
      if (j0 == 0) j0 <- lambda.length1
      k0 <- (location-j0) / lambda.length1 + 1
      lambda <- exp(c(Lambda1[j0], Lambda1[j0]))
    } ## end of search.grid

    if (search.grid == F) {
      fit <- stats::optim(
        0,
        fbps_gcv,
        method = method,
        control = control,
        lower = lower[1],
        upper = upper[1]
      )

      fit <- stats::optim(
        c(fit$par, fit$par),
        fbps_gcv,
        method = method,
        control = control,
        lower = lower[1],
        upper = upper[1]
      )
      if (fit$convergence>0) {
        expression <- paste(
          "Smoothing failed! The code is:",
          fit$convergence
        )
        print(expression)
      }
      lambda <- exp(fit$par)
    } ## end of optim

  } ## end of finding smoothing parameters
  lambda = rep(lambda, 2)[1:2]

  return(fbps_est(log(lambda)))
}
