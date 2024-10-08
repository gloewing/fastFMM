#' Generic "G_generate" G matrix setup
#'
#' Creates the design matrix that allows for the estimation of the G matrix.
#' Output related to "G_estimate" dispatch.
#'
#' @param mum An object that inherits from the "massmm" class
#' @param ... Additional arguments
#'
#' @return A design matrix for G estimation.

G_generate <- function(mum, ...) {
  UseMethod("G_generate")
}

#' Generate the G matrix in non-concurrent cases
#'
#' Creates a design matrix from the massive univariate step, assuming the model
#' is not concurrent.
#'
#' @param mum The massively univariate model output. Contains relevant data,
#' such as the transposed Ztlist from the univariate fit.
#' @param MoM indicator for type of MoM estimation
#'
#' @return A list with the new Z, the original Ztlist, and indices

G_generate.massmm <- function(mum, MoM) {

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
#' @param mum The massively univariate model output with class "massmm_conc".
#' Contains relevant data, such as the transposed Ztlist from the univariate fit.
#' @param MoM indicator for type of MoM estimation
#' @param i Integer for first index of the covariance
#' @param j Integer for second index of the covariance
#'
#' @return A list with the new Z, the original Ztlist, and indices
#'
#' @importFrom Matrix rowSums

# AX: Add back all_crossterms
G_generate.massmm_conc <- function(mum, i, j) {

  # Z_lst is the ZTlist (transposed) output from:
  #   sapply(getME(fit_uni, "Ztlist"), function(x) t(x) )
  # RE_table is a table from VarCorr(X) where X is a lme4 "lmerMod" class object
  # ID is the name of the ID factor (which determines what we can sum across)
  # assumes the names of HHat are generated from this same table in same order

  RE_table <- mum$varcorr_df
  ID <- mum$subj_id

  # 1) Gather the appropriate Z matrices

  # AX: It's possible the base rowSums will also work
  Z_list_i <- mum$ztlist[[i]]
  Z_orig_i <- sapply(Z_list_i, Matrix::rowSums)

  Z_list_j <- mum$ztlist[[j]]
  Z_orig_j <- sapply(Z_list_j, Matrix::rowSums)

  #1) concatenate Z_orig
  # list where each element will be the design submatrix for corresponding term
  Z_orig <- Z_lst
  z_names <- names(Z_lst)

  Z <- vector(length = nrow(RE_table), "list")

  for (k in 1:nrow(RE_table)) {
    # iterate through covariance term names (i.e., random effect terms)

    cross_term <- !is.na(RE_table$var2)[k]  # cross term (covariance term)
    re_name <-  RE_table[k, 1] # random effects

    # check if interaction
    if (grepl(":", re_name, fixed = TRUE)) {
      re_interact <- TRUE # interaction of random effects
      ID_flag <- FALSE # this is always false for interactions of random effects
    } else {
      re_interact <- FALSE # interaction of random effects
      # this determines whether the main subject/ID variable is triggered;
      # indicates we should rowSum across all columns associated with ID
      ID_flag <- ifelse(re_name == ID, TRUE, FALSE)
    }

    # either a blank or the name of the last variable
    v2 <- ifelse(is.na(RE_table$var2[k]), "", paste0("_", RE_table$var2[k]))
    # intercept term so does not require squaring
    intrcpt <- ifelse(RE_table$var1[k] == "(Intercept)", TRUE, FALSE )

    # check to see if this is the cross-term between two random effects (e.g., intercept x slope)
    if (!cross_term) {
      # NOT a cross-term (covariance term)

      # find grp -- var1 combination that matches Z_lst names (z_names)
      var1 <- RE_table$var1[k]
      nm <- paste0(re_name, ".",  var1) # name of Z_lst (outputted by lme4 getME() )
      zlst_idx <- which(z_names == nm) # find index of zlst element that has appropriate design matrix

      Z[[k]] <- Z_list_i[[zlst_idx]] * Z_list_j[[zlst_idx]]

    } else {
      # Cross term

      ## since cross term is element-wise product between a random intercept and a random slope, the entries are only non-zero
      ## when the random intercept is 1 (which are the same for the same factor RE), so we just need to multiple
      ## the random slope by 2 (to emulate the cross term), all other terms will be 0 (so we can avoid those)

      # find grp -- var1 combination that matches Z_lst names (z_names)
      rand_slope <- RE_table$var2[k] # the cross terms do not use var1, they only use var2 for znames
      nm <- paste0( re_name, ".",  rand_slope ) # name of Z_lst (outputted by lme4 getME() )
      rand_slope_idx <- which(z_names == nm) # find index of matrix that has appropriate design matrix

      # TODO: Check if this is equivalent to Z_lst[[rand_slope_idx]] * Z_lst [[intrcpt]] * 2
      Z[[k]] <- Z_list_i[[rand_slope_idx]] + Z_list_j[[rand_slope_idx]]
    }

    # ID flag -- if main ID variable is the only random effect factor for row j of RE_table (like (1 | ID  )   or (variable | ID), then these submatrices are summed across columns )
    # sum across columns

    if (ID_flag) {
      Z[[k]] <- matrix(Matrix::rowSums(Z[[k]]), ncol = 1)
      colnames(Z[[k]]) <- paste0(RE_table$grp[k], "_", RE_table$var1[k], v2) # name column
      idx_vec[k] <- 1
    } else {
      idx_vec[k] <- ncol(Z[[k]])         # number of columns in matrix
      colnames(Z[[k]]) <- paste0(RE_table$grp[k], "_", RE_table$var1[k], v2, "_", 1:ncol(Z[[k]]))   # name columns
    }
  }

  Z <- do.call(cbind, Z) # concatenate

  # TODO: Add cross-terms for multiple covariates/columns same random effects
  # Probably need to doctor the data a bit to test this
  # AX: Add all_crossterms instead of addition

  return(
    list(
      Z = Z,
      Z_orig = list(Z_orig_i, Z_orig_j),
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
