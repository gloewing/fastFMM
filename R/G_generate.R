#' Creates the design matrix that allows for estimation of G
#'
#' The function `G_estimate` uses a MoM method,
#' and `G_estimate_randint` is a special case of `G_estimate`.
#'
#' Because `G_estimate`
#'
#' @param data Data frame that contains the predictors and outcome
#' @param Z_lst Transposed list of Z matrices from the univariate fits
#' @param RE_table Table of random effects and interactions, generated from the
#' `lmerMod` object
#' @param ID Name of the ID factor, assuming names of `HHat` are generated from
#' the same table in the same order
#'
#' @return List containing Z matrices and indices (unsure)
#'
#' @import lme4
#'
#' @export

G_generate <- function(data, Z_lst, RE_table, ID ="id"){
  # data fed to fui function
  # Z_lst is the ZTlist (transposed) output from:
  #   sapply(getME(fit_uni, "Ztlist"), function(x) t(x) )
  # RE_table is a table from VarCorr(X) where X is a lme4 "lmerMod" class object
  # ID is the name of the ID factor (which determines what we can sum across)
  # assumes the names of HHat are generated from this same table in same order

  #1) concatenate Z_orig
  Z_orig <- Z_lst # list where each element will be design submatrix for corresponding term
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
  Z <- vector(length = nrow(RE_table), "list") # list where each element will be design submatrix for corresponding term
  idx_vec <- vector(length = nrow(RE_table)) # vector where each element are indices of eventual Z matrix corresponding to each term in HHat (for summing in OLS)

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
      ID_flag <- ifelse (re_name == ID, TRUE, FALSE) # this determines whether the main subject/ID variable is triggered -- indicates we should rowSum across all columns associated with ID
    }

    v2 <- ifelse (is.na(RE_table$var2[k]), "", paste0("_", RE_table$var2[k])) # either a blank or the name of the last variable
    intrcpt <- ifelse ( RE_table$var1[k] == "(Intercept)", TRUE, FALSE )  # intercept term so does not require squaring

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
  Z_orig <- sapply(Z_orig, rowSums) # sum all sub matrices for random effects design matrix
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
