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

G_generate <- function(data, Z_lst, RE_table, ID ="id"){
  # data fed to fui function
  # Z_lst is the ZTlist (transposed) output from: sapply(getME(fit_uni, "Ztlist"), function(x) t(x) )
  # RE_table is a table from VarCorr( X ) where X is a lme4 "lmerMod" class object
  # ID is the name of the ID factor (which determines what we can sum across)
  # assumes the names of HHat are generated from this same table in same order

  #1) concatenate Z_orig
  Z_orig <- Z_lst # list where each element will be design submatrix for corresponding term
  z_names <- names(Z_lst)

  # sum across variables -- this comes from derivations in Overleaf for method of moments estimator of G(s_1, s_2)
  Z_orig <- sapply(Z_orig, rowSums) # sum all sub matrices for random effects design matrix
  colnames(Z_orig) <- z_names

  #2) prepare design matrix for regression
  Z <- vector(length = nrow(RE_table), "list") # list where each element will be design submatrix for corresponding term
  idx_vec <- vector(length = nrow(RE_table)) # vector where each element are indices of eventual Z matrix corresponding to each term in HHat (for summing in OLS)

  for(j in 1:nrow(RE_table)){
    # iterate through covariance term names (i.e., random effect terms)

    cross_term <- !is.na(RE_table$var2)[j]  # cross term (covariance term)
    re_name <-  RE_table[j, 1] # random effects

    # check if interaction
    if(grepl(":", re_name, fixed = TRUE)){
      re_interact <- TRUE # interaction of random effects
      ID_flag <- FALSE # this is always false for interactions of random effects
    }else{
      re_interact <- FALSE # interaction of random effects
      ID_flag <- ifelse(re_name == ID, TRUE, FALSE) # this determines whether the main subject/ID variable is triggered -- indicates we should rowSum across all columns associated with ID
    }

    v2 <- ifelse(is.na(RE_table$var2[j]), "", paste0("_", RE_table$var2[j])) # either a blank or the name of the last variable
    intrcpt <- ifelse( RE_table$var1[j] == "(Intercept)", TRUE, FALSE )  # intercept term so does not require squaring

    # check to see if this is the cross-term between two random effects (e.g., intercept x slope)
    if(!cross_term ){
      # NOT a cross-term (covariance term)

      # find grp -- var1 combination that matches Z_lst names (z_names)
      var1 <- RE_table$var1[j]
      nm <- paste0( re_name, ".",  var1 ) # name of Z_lst (outputted by lme4 getME() )
      zlst_idx <- which(z_names == nm) # find index of zlst element that has appropriate design matrix

      if(intrcpt){
        # intercept term (i.e., does not require squaring elements) since indicators squared are just indicators
        Z[[j]] <- Z_lst[[zlst_idx]]
      }else{
        # not an intercept term (corresponds to random slope so requires squaring elements) -- see below for why we can square instead of doing actual element-wise produt (all other product terms are zero-ed out)
        Z[[j]] <- (Z_lst[[zlst_idx]])^2
      }

    } else {
      # cross term

      ## since cross term is element-wise product between a random intercept and a random slope, the entries are only non-zero
      ## when the random intercept is 1 (which are the same for the same factor RE), so we just need to multiple
      ## the random slope by 2 (to emulate the cross term), all other terms will be 0 (so we can avoid those)

      # find grp -- var1 combination that matches Z_lst names (z_names)
      rand_slope <- RE_table$var2[j] # the cross terms do not use var1, they only use var2 for znames
      nm <- paste0( re_name, ".",  rand_slope ) # name of Z_lst (outputted by lme4 getME() )
      rand_slope_idx <- which(z_names == nm) # find index of matrix that has appropriate design matrix

      Z[[j]] <- Z_lst[[ rand_slope_idx ]] * 2 # element of Z_lst corresponding to appropriate random slope (scale by 2 to emulate cross term)
    }

    # ID flag -- if main ID variable is the only random effect factor for row j of RE_table (like (1 | ID  )   or (variable | ID), then these submatrices are summed across columns )
    # sum across columns

    if (ID_flag) {
      Z[[j]] <- matrix( rowSums(Z[[j]]), ncol = 1)
      colnames(Z[[j]]) <- paste0(
        RE_table$grp[j], "_",
        RE_table$var1[j],
        v2
      ) # name column
      idx_vec[j] <- 1
    } else {
      idx_vec[j] <- ncol(Z[[j]])         # number of columns in matrix
      colnames(Z[[j]]) <- paste0(
        RE_table$grp[j], "_",
        RE_table$var1[j], v2, "_",
        1:ncol(Z[[j]])
      )   # name columns
    }

  }

  idx_vec <- c(0, cumsum(idx_vec))
  idx_lst <- sapply(
    seq_along(1:nrow(RE_table)),
    function(x) (idx_vec[x] + 1):(idx_vec[x+1])
  ) # column indices
  Z <- do.call(cbind, Z) # concatenate

  return(
    list(
      Z = Z,
      Z_orig = Z_orig,
      idx_lst = idx_lst,
      idx_vec = idx_vec[-1]
    )
  )

}
