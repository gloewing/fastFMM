#' Generic "G_generate" G matrix setup
#'
#' Creates the design matrix that allows for the estimation of the G matrix.
#' Output related to "G_estimate" dispatch.
#'
#' @param fmm "fastFMM" object
#' @param ... Additional arguments
#'
#' @return A design matrix for G estimation.
#' @export

G_generate <- function(fmm, ...) {
  UseMethod("G_generate")
}

#' Generate the G matrix in non-concurrent cases
#'
#' Creates a design matrix from the massive univariate step, assuming the model
#' is not concurrent.
#'
#' @param fmm "fastFMM" object. Ignored except for dispatch.
#' @param mum The massively univariate model output. Contains relevant data,
#' such as the transposed Ztlist from the univariate fit.
#' @param MoM indicator for type of MoM estimation
#'
#' @method G_generate fastFMM
#' @return A list with the new Z, the original Ztlist, and indices
#' @export

G_generate.fastFMM <- function(fmm, mum, MoM) {

  # Z_lst is the ZTlist (transposed) output from:
  #   sapply(getME(fit_uni, "Ztlist"), function(x) t(x) )
  # RE_table is a table from VarCorr(X) where X is a lme4 "lmerMod" class object
  # ID is the name of the ID factor (which determines what we can sum across)
  # assumes the names of HHat are generated from this same table in same order

  Z_lst <- mum$ztlist
  RE_table <- mum$varcorr_df
  ID <- mum$subj_id

  #1) concatenate Z_orig
  # list where each element will be design submatrix for corresponding term
  Z_orig <- Z_lst
  z_names <- names(Z_lst)

  # Use summed rows if MoM == 1; o/w, use unsummed version
  if (MoM == 1) {
    Z_list_i <- Z_list_j <- lapply(
      Z_lst, function(x) matrix(rowSums(x), ncol = 1)
    )
  } else {
    Z_list_i <- Z_list_j <- Z_lst
  }


  #2) prepare design matrix for regression
  # list where each element will be design submatrix for corresponding term
  Z <- vector(length = nrow(RE_table), "list")
  # vector where each element are indices of eventual Z matrix corresponding to
  # each term in HHat (for summing in OLS)
  idx_vec <- vector(length = nrow(RE_table))

  # iterate through covariance term names (i.e., random effect terms)
  for (k in 1:nrow(RE_table)) {
    nm1 <- paste0(RE_table[k, 1], '.', RE_table$var1[k])
    idx1 <- which(z_names == nm1)

    # When MoM == 1, avoid using all_crossterms and take the product
    if (MoM == 1) {
      if (is.na(RE_table$var2)[k]) { # not cross-term
        Z[[k]] <- matrix(
          Z_list_i[[idx1]][, 1] * Z_list_j[[idx1]][, 1],
          ncol = 1
        )
      } else { # is cross-term
        nm2 <- paste0(RE_table[k, 1], '.', RE_table$var2[k])
        idx2 <- which(z_names == nm2)
        Z[[k]] <- matrix(
          Z_list_i[[idx1]][, 1] * Z_list_j[[idx2]][, 1] * 2,
          ncol = 1
        )
      }
    } else { # MoM == 2
      if (is.na(RE_table$var2)[k]) { # not crossterm
        Z[[k]] <- all_crossterms(Z_list_i[[idx1]], Z_list_j[[idx1]])
      } else { # is crossterm
        nm2 <- paste0(RE_table[k, 1], '.', RE_table$var2[k])
        idx2 <- which(z_names == nm2)
        Z[[k]] <- all_crossterms(Z_list_i[[idx1]], Z_list_j[[idx2]]) * 2
      }
    }

    re_name <-  RE_table[k, 1] # random effects

    # check if interaction
    if (grepl(":", re_name, fixed = TRUE)) {
      re_interact <- TRUE # interaction of random effects
      ID_flag <- FALSE # this is always false for interactions of random effects
    } else {
      re_interact <- FALSE # interaction of random effects
      # This determines whether the main subject/ID variable is triggered:
      # Indicates we should rowSum across all columns associated with ID
      ID_flag <- ifelse (re_name == ID, TRUE, FALSE)
    }

    # either a blank or the name of the last variable
    v2 <- ifelse (is.na(RE_table$var2[k]), "", paste0("_", RE_table$var2[k]))
    # intercept term so does not require squaring
    intrcpt <- ifelse ( RE_table$var1[k] == "(Intercept)", TRUE, FALSE )

    colnames(Z[[k]]) <- paste0(
      RE_table$grp[k], "_", RE_table$var1[k], v2, "_", 1:ncol(Z[[k]])
    )

    if (ID_flag & MoM == 2)
      idx_vec[k] <- 1
    else
      idx_vec[k] <- ncol(Z[[k]])

  }

  idx_vec <- c(0, cumsum(idx_vec))
  # column indices
  idx_lst <- sapply(
    seq_along(1:nrow(RE_table)),
    function(x) (idx_vec[x] + 1):(idx_vec[x+1])
  )
  Z <- do.call(cbind, Z)

  # sum across variables
  # this comes from derivations in Overleaf for method of moments estimator of G(s_1, s_2)
  # sum all sub matrices for random effects design matrix
  Z_orig <- sapply(Z_orig, rowSums)
  colnames(Z_orig) <- z_names

  return(
    list(
      Z = Z,
      Z_orig = Z_orig,
      idx_lst = idx_lst,
      idx_vec = idx_vec[-1]
    )
  )
}

#' Generate the G matrix in concurrent cases
#'
#' Creates a design matrix from the massive univariate step, assuming the model
#' is concurrent. Unlike the non-concurrent case, there is curently no encoding
#' for different MoM estimators.
#'
#' @param fmm "fastFMM" object. Ignored except for dispatch.
#' @param mum The massively univariate model output. Contains relevant data,
#' such as the transposed Ztlist from the univariate fit.
#' @param MoM indicator for type of MoM estimation
#' @param i Integer for first index of the covariance
#' @param j Integer for second index of the covariance
#'
#' @method G_generate fastFMMconc
#' @return A list with the new Z, the original Ztlist, and indices
#'

#' @importFrom Matrix rowSums
#' @export

# AX: Add back all_crossterms
G_generate.fastFMMconc <- function(fmm, mum, i, j, MoM = 1) {

  # Z_lst is the ZTlist (transposed) output from:
  #   sapply(getME(fit_uni, "Ztlist"), function(x) t(x) )
  # RE_table is a table from VarCorr(X) where X is a lme4 "lmerMod" class object
  # ID is the name of the ID factor (which determines what we can sum across)
  # assumes the names of HHat are generated from this same table in same order

  # Choice of i or j for the RE table is arbitrary
  RE_table <- mum$varcorr_df[[i]]
  ID <- mum$subj_id

  # 1) Gather the appropriate Z matrices

  # Note these are rowSums
  Z_list_i <- mum$ztlist[[i]]
  Z_list_j <- mum$ztlist[[j]]

  # list where each element will be the design submatrix for corresponding term
  z_names_i <- colnames(Z_list_i)
  z_names_j <- colnames(Z_list_j)
  # RE_table rows should be consistent, so i or j is arbitrary
  Z <- vector(length = nrow(RE_table), "list")
  idx_vec <- vector(length = nrow(RE_table))

  for (k in 1:nrow(RE_table)) {
    re_name <-  RE_table[k, 1] # random effects
    var1_i <- RE_table$var1[k]
    var1_j <- gsub(paste0("\\_", i, "$"), paste0("\\_", j), var1)
    nm1_i <- paste0(re_name, ".",  var1_i) # name of Z_lst (outputted by lme4 getME() )
    nm1_j <- paste0(re_name, ".", var1_j)

    if (is.na(RE_table$var2)[k]) { # not cross-term
      Z[[k]] <- matrix(
        Z_list_i[, nm1_i] * Z_list_j[, nm1_j],
        ncol = 1
      )
    } else { # is cross-term
      # Pick nm2 from the j index
      var2_j <- gsub(paste0("\\_", i, "$"), paste0("\\_", j), RE_table$var2[k])
      nm2_j <- paste0(RE_table[k, 1], '.', var2_j)
      Z[[k]] <- matrix(
        Z_list_i[, nm1_i] * Z_list_j[, nm2_j] * 2,
        ncol = 1
      )
    }
    idx_vec[k] <- ncol(Z[[k]])
  }

  idx_vec <- c(0, cumsum(idx_vec))
  # column indices
  idx_lst <- sapply(
    seq_along(1:nrow(RE_table)),
    function(x) (idx_vec[x] + 1):(idx_vec[x+1])
  )
  Z <- do.call(cbind, Z)

  return(
    list(
      Z = Z,
      # Z_orig is already gone because of rowSums previously
      Z_orig = NULL,
      ztlist = list(Z_list_i, Z_list_j),
      idx_lst = idx_lst,
      idx_vec = idx_vec[-1],
      RE_table = RE_table
    )
  )
}

#' Create crossterms from two matrices
#'
#' A helper function for `G_generate` that produces cross-terms.
#'
#' @param Z_i Matrix
#' @param Z_j Matrix
#' @param make_sparse Boolean for whether to output a sparse matrix.
#' Default is `TRUE`.
#'
#' @return Matrix of cross-terms between `Z_i` and `Z_j`.
#'
#' @import Matrix

all_crossterms <- function(Z_i, Z_j, make_sparse = TRUE) {

  if (!identical(dim(Z_i), dim(Z_j)))
    stop('In G_generate, matrices must have the same dimension for cross-term.')

  n <- ncol(Z_i)
  # Get combinations without repeats
  temp <- cbind(rep(1:n, n:1), do.call(c, sapply(1:n, function(x) x:n)))
  xterm <- sapply(
    1:nrow(temp),
    function(idx) Z_i[, temp[idx, 1]] * Z_j[, temp[idx, 2]]
  )

  if (make_sparse)
    return(Matrix(xterm, sparse = T))
  else
    return(xterm)
}
