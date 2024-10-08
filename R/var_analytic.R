#' Analytic variance calculation
#'
#' @param mum Massively univariate model output of class "massmm"
#' @param betaHat Matrix of smoothed coefficients
#' @param data Data frame of values to fit
#' @param L integer number of points on functional domain
#' @param MoM integer MoM estimator
#' @param non_neg Integer indicator of non-negativity calculation
#' @param nknots_cov Integer number of knots for splines
#' @param parallel Boolean indicator for parallel processing
#' @param silent Logical, suppresses messages when `TRUE`. Passed from `fui`.
#'
#' @return List of final outputs of `fui`
#'
#' @importFrom stats smooth.spline quantile
#' @importFrom Rfast spdinv rowMaxs
#' @importFrom Matrix tcrossprod crossprod Diagonal
#' @importFrom parallel mclapply
#' @importFrom methods new
#' @importFrom mvtnorm rmvnorm

var_analytic <- function(
  mum,
  betaHat,
  data,
  L,
  MoM,
  non_neg,
  nknots_cov,
  parallel,
  silent
) {

# 0 Preparation of H, R #####################################################

  # Variance estimates of random components: H(s), R(s)
  if (!silent) print("Step 3.0: Preparation")

  randintercept <- mum$randintercept
  argvals <- mum$argvals
  p <- nrow(mum$betaTilde)

  HHat <- t(
    apply(
      mum$sigmausqHat, 1,
      function(b) stats::smooth.spline(x = mum$argvals, y = b)$y
    )
  )
  ind_var <- which(grepl("var", rownames(HHat)) == TRUE)
  HHat[ind_var, ][which(HHat[ind_var, ] < 0)] <- 0
  RHat <- t(
    apply(
      mum$sigmaesqHat,
      1,
      function(b) stats::smooth.spline(x = mum$argvals, y = b)$y
    )
  )
  RHat[which(RHat < 0)] <- 0

  # 1 MoM estimator of G #######################################################

  # 1.1 Preparation A: Random intercept flag ===================================

  # Includes potential NNLS correction for diagonals (for variance terms)

  if (randintercept) {

    GTilde <- G_estimate_randint(
      data = data,
      L = L,
      out_index = mum$out_index,
      designmat = mum$designmat,
      betaHat = betaHat,
      silent = silent
    )

    if (!silent) {
      print("Step 3.1: MoM estimator of G")
      print("Step 3.1.1: Preparation A")
    }

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

    # 1.1 Preparation B ========================================================

    ## if more random effects than random intercept only
    if (!silent)
      print("Step 3.1.1: Preparation B")

    # generate design matrix for G(s_1, s_2)
    # Method of Moments Linear Regression calculation

    ## Method of Moments estimator for G() with potential NNLS correction
    # for diagonals (for variance terms)

    G_est_out <- G_estimate(
      mum = mum,
      data = data,
      L = L,
      betaHat = betaHat,
      HHat = HHat,
      non_neg = non_neg,
      MoM = MoM,
      silent = silent
    )
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

  ## Obtain the corresponding rows of each subject
  obs.ind <- list()
  group <- mum$subj_id
  ID.number <- unique(data[, group])
  for(id in ID.number) {
    obs.ind[[as.character(id)]] <- which(data[, group] == id)
  }

  ## Concatenate vector of 1s to Z because used that way below
  HHat_trim <- NA
  if (!randintercept) {
    Z <- mum$ztlist
    qq <- ncol(Z)
    HHat_trim <- array(NA, c(qq, qq, L)) # array for Hhat
  }

  ## Create var.beta.tilde.theo to store variance of betaTilde
  var.beta.tilde.theo <- array(NA, dim = c(p, p, L))
  ## Create XTVinvZ_all to store all XTVinvZ used in the
  # covariance calculation
  XTVinvZ_all <- vector(length = L, "list")
  ## arbitrarily start find indices
  resStart <- cov_organize_start(HHat[,1])
  res_template <- resStart$v_list_template # index template
  template_cols <- ncol(res_template)

  # AX: Take this out of the function
  ## Calculate Var(betaTilde) for each location
  parallel_fn <- function(s){
    V.subj.inv <- c()
    ## Invert each block matrix, then combine them
    if (!randintercept) {
      cov.trimmed <- eigenval_trim(
        matrix(c(HHat[, s], 0)[res_template], template_cols)
      )
      HHat_trim[, , s] <- cov.trimmed
    }

    # store XT * Vinv * X
    XTVinvX <- matrix(0, nrow = p, ncol = p)
    # store XT * Vinv * Z
    XTVinvZ_i <- vector(length = length(ID.number), "list")

    ## iterate for each subject
    for (id in ID.number) {
      subj_ind <- obs.ind[[as.character(id)]]
      subj.ind <- subj_ind

      if (randintercept) {
        Ji <- length(subj_ind)
        V.subj <- matrix(HHat[1, s], nrow = Ji, ncol = Ji) +
          diag(RHat[s], Ji)
      } else {
        V.subj <- Z[subj_ind, , drop = FALSE] %*%
          tcrossprod(cov.trimmed, Z[subj_ind, , drop = FALSE]) +
          diag(RHat[s], length(subj_ind))
      }

      # Rfast provides a faster matrix inversion
      V.subj.inv <- as.matrix(Rfast::spdinv(V.subj))

      XTVinvX <- XTVinvX +
        crossprod(matrix(designmat[subj.ind,], ncol = p), V.subj.inv) %*%
        matrix(designmat[subj.ind,], ncol = p)

      if (randintercept) {
        XTVinvZ_i[[as.character(id)]] <- crossprod(
          matrix(designmat[subj.ind,], ncol = p),
          V.subj.inv
        ) %*% matrix(1, nrow = Ji, ncol = 1)
      } else {
        XTVinvZ_i[[as.character(id)]] <- crossprod(
          matrix(designmat[subj.ind,], ncol = p),
          V.subj.inv
        ) %*% Z[subj.ind,,drop = FALSE]
      }
    }

    var.beta.tilde.theo[, , s] <- as.matrix(Rfast::spdinv(XTVinvX))

    return(
      list(
        XTViv = XTVinvZ_i,
        var.beta.tilde.theo = var.beta.tilde.theo[,,s]
      )
    )
  }

  # Fit massively univariate mixed models
  if (parallel) {
    # Windows is incompatible with mclapply
    if(.Platform$OS.type == "windows") {
      # Provide necessary objects to the cluster
      clusterExport(
        cl = cl,
        varlist = c(
          "eigenval_trim",
          "HHat",
          "s",
          "res_template",
          "template_cols",
          "p",
          "ID.number",
          "obs.ind",
          "id",
          "RHat"
        ),
        # Look within the function's environment
        envir = environment()
      )
      if (!randintercept)
        clusterExport(cl, "Z", envir = environment())

      # tryCatch ensures the cluster is stopped even during errors
      massVar <- tryCatch(
        parLapply(
          cl = cl,
          X = argvals,
          fun = parallel_fn
        ),
        warning = function(w) {
          print(
            paste0(
              "Warning when parallelizing variance calculation:", "\n",
              w
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
        FUN = parallel_fn,
        mc.cores = num_cores
      )
    }
  } else {
    massVar <- lapply(
      argvals,
      parallel_fn
    )
  }

  XTVinvZ_all <- lapply(argvals, function(s)
    massVar[[s]]$XTViv[lengths(massVar[[s]]$XTViv) != 0])
  var.beta.tilde.theo <- lapply(argvals, function(s)
    massVar[[s]]$var.beta.tilde.theo)
  var.beta.tilde.theo <- simplify2array(var.beta.tilde.theo) %>%
    array(dim = c(p, p, L))

  suppressWarnings(rm(massVar, resStart, res_template, template_cols))

  # Calculate the inter-location covariance of betaTilde:
  # Cov(betaTilde(s_1), betaTilde(s_2))
  ## Create cov.beta.tilde.theo to store covariance of betaTilde
  cov.beta.tilde.theo <- array(NA, dim = c(p,p,L,L))
  if (randintercept) {
    resStart <- cov_organize_start(GHat[1,2]) # arbitrarily start
  } else {
    resStart <- cov_organize_start(GHat[,1,2]) # arbitrarily start
  }

  # 2.1 First step =============================================================

  if (!silent) print("Step 3.2.1: First step")

  res_template <- resStart$v_list_template # index template
  template_cols <- ncol(res_template)
  ## Calculate Cov(betaTilde) for each pair of location
  for (i in 1:L) {
    for (j in i:L) {
      V.cov.subj <- list()
      tmp <- matrix(0, nrow = p, ncol = p) ## store intermediate part
      if (randintercept) {
        G_use <- GHat[i, j]
      } else {
        G_use <- eigenval_trim(
          matrix(c(GHat[, i, j], 0)[res_template], template_cols)
        )
      }

      for (id in ID.number) {
        tmp <- tmp +
          XTVinvZ_all[[i]][[as.character(id)]] %*%
          tcrossprod(G_use, XTVinvZ_all[[j]][[as.character(id)]])
      }

      ## Calculate covariance using XTVinvX and tmp to save memory
      cov.beta.tilde.theo[, , i, j] <- var.beta.tilde.theo[, , i] %*%
        tmp %*%
        var.beta.tilde.theo[, , j]
    }
  }
  suppressWarnings(
    rm(
      V.subj,
      V.cov.subj,
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
  var.beta.tilde.s <- array(NA, dim = c(L, L, p))
  for(j in 1:p) {
    for(r in 1:L) {
      for(t in 1:L) {
        if (t == r) {
          var.beta.tilde.s[r,t,j] <- var.beta.tilde.theo[j, j, r]
        } else {
          var.beta.tilde.s[r,t,j] <- cov.beta.tilde.theo[
            j, j, min(r, t), max(r, t)
          ]
        }
      }
    }
  }

  # Calculate the inter-location covariance of betaHat:
  # Cov(betaHat(s_1), betaHat(s_2))
  var.beta.hat <- array(NA, dim = c(L, L, p))
  for(r in 1:p) {
    M <- B %*%
      Matrix::tcrossprod(solve(Matrix::crossprod(B) + lambda[r] * S), B) +
      matrix(1/L, nrow = L, ncol = L)
    var.raw <- M %*% Matrix::tcrossprod(var.beta.tilde.s[, , r], M)
    ## trim eigenvalues to make final variance matrix PSD
    var.beta.hat[,,r] <- eigenval_trim(var.raw)
  }
  betaHat.var <- var.beta.hat ## final estimate

  # Obtain qn to construct joint CI
  qn <- rep(0, length = nrow(betaHat))
  N <- 10000 ## sample size in simulation-based approach
  zero_vec <- rep(0, length(betaHat[1,]))
  set.seed(seed)

  for(i in 1:length(qn)) {
    Sigma <- betaHat.var[, , i]
    sqrt_Sigma <- sqrt(diag(Sigma))
    S_scl <- Matrix::Diagonal(x = 1 / sqrt_Sigma)
    Sigma <- as.matrix(S_scl %*% Sigma %*% S_scl)
    # x_sample <- abs(FastGP::rcpp_rmvnorm_stable(N, Sigma, zero_vec))
    x_sample <- abs(mvtnorm::rmvnorm(N, zero_vec, Sigma))
    un <- Rfast::rowMaxs(x_sample, value = TRUE)
    qn[i] <- stats::quantile(un, 0.95)
  }

  # Decide whether to return design matrix or just set it to NULL
  if (!design_mat) designmat <- NA
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
      betaHat.var = betaHat.var,
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
      argvals = mum$argvals,
      randeffs = mum$randeffs,
      se_mat = mum$se_mat
    )
  )

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

eigenval_trim <- function(V) {
  ## trim non-positive eigenvalues to ensure positive semidefinite
  edcomp <- eigen(V, symmetric = TRUE)
  eigen.positive <- which(edcomp$values > 0)
  q <- ncol(V)

  if (length(eigen.positive) == q) {
    # nothing needed here because matrix is already PSD
    return(V)
  } else if (length(eigen.positive) == 0) {
    return(tcrossprod(edcomp$vectors[, 1]) * edcomp$values[1])
  } else if (length(eigen.positive) == 1) {
    return(
      tcrossprod(as.vector(edcomp$vectors[,1])) *
        as.numeric(edcomp$values[1])
    )
  } else {
    # sum of outer products of eigenvectors, scaled by eigenvalues for all
    # positive eigenvalues
    return(
      matrix(
        edcomp$vectors[, eigen.positive] %*%
          tcrossprod(
            diag(edcomp$values[eigen.positive]),
            edcomp$vectors[,eigen.positive]
          ),
        ncol = q
      )
    )
  }
}
