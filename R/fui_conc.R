#' Concurrent FUI
#'
#' The concurrent version of `fui.R`.
#'
#' Fit a function-on-scalar regression model for longitudinal
#' functional outcomes and scalar predictors using the Fast Univariate
#' Inference (FUI) approach (Cui et al. 2022).
#'
#' The FUI approach comprises of three steps:
#' \enumerate{
#' \item Fit a univariate mixed model at each location of the functional domain,
#' and obtain raw estimates from massive models;
#' \item Smooth the raw estimates along the functional domain;
#' \item Obtain the pointwise and joint confidence bands using an analytic
#' approach for Gaussian data or Bootstrap for general distributions.
#' }
#'
#' For more information on each step, please refer to the FUI paper
#' by Cui et al. (2022).
#'
#' @param formula Two-sided formula object in lme4 formula syntax.
#' The difference is that the response need to be specified as a matrix
#' instead of a vector. Each column of the matrix represents one location
#' of the longitudinal functional observations on the domain.
#' @param data A data frame containing all variables in formula
#' @param family GLM family of the response. Defaults to \code{gaussian}.
#' @param var Logical, indicating whether to calculate and return variance
#' of the coefficient estimates. Defaults to \code{TRUE}.
#' @param analytic Logical, indicating whether to use the analytic inference
#' approach or bootstrap. Defaults to \code{TRUE}.
#' @param parallel Logical, indicating whether to do parallel computing.
#' Defaults to \code{FALSE}.
#' @param silent Logical, indicating whether to show descriptions of each step.
#' Defaults to \code{FALSE}.
#' @param argvals A vector containing locations of observations on the
#' functional domain. If not specified, a regular grid across the range of
#' the domain is assumed. Currently only supported for bootstrap (\code{analytic=FALSE}).
#' @param nknots_min Minimal number of knots in the penalized smoothing for the regression coefficients.
#' Defaults to \code{NULL}, which then uses L/2 where L is the dimension of the functional domain.
#' @param nknots_min_cov Minimal number of knots in the penalized smoothing for the covariance matrices.
#' Defaults to \code{35}.
#' @param smooth_method How to select smoothing parameter in step 2. Defaults to \code{"GCV.Cp"}
#' @param splines Spline type used for penalized splines smoothing. We use the same syntax
#' as the mgcv package. Defaults to \code{"tp"}
#' @param design_mat Logical, indicating whether to return the design matrix.
#' Defaults to \code{FALSE}
#' @param residuals Logical, indicating whether to save residuals from
#' unsmoothed LME. Defaults to \code{FALSE}.
#' @param G_return Logical, indicating whether to return (smoothed and trimmed)
#' G = Cov(u(s_t), u(s_l)). Defaults to \code{FALSE}.
#' @param num_boots Number of samples when using bootstrap inference. Defaults to 500.
#' @param boot_type Bootstrap type (character): "cluster", "case", "wild", "reb",
#' "residual", "parametric", "semiparametric". \code{NULL} defaults to "cluster"
#' for non-gaussian responses and "wild" for gaussian responses. For small cluster
#' (n<=10) gaussian responses, defaults to "reb"
#' @param seed Numeric value used to make sure bootstrap replicate (draws) are
#' correlated across functional domains for certain bootstrap approach
#' @param subj_ID Name of the variable that contains subject ID.
#' @param num_cores Number of cores for parallelization. Defaults to 1.
#' @param caic Logical, indicating whether to calculate cAIC. Defaults to \code{FALSE}.
#' @param REs Logical, indicating whether to return random effect estimates.
#' Defaults to \code{FALSE}.
#' @param non_neg 0 - no non-negativity constrains, 1 - non-negativity
#' constraints on every coefficient for variance, 2 - non-negativity on
#' average of coefficents for 1 variance term. Defaults to 0.
#' @param MoM Method of moments estimator. Default to 2. 1 should only be used
#' for extremely large datasets.
#'
#' @return A list containing:
#' \item{betaHat}{Estimated functional fixed effects}
#' \item{argvals}{Location of the observations}
#' \item{betaHat.var}{Variance estimates of the functional fixed effects (if specified)}
#' \item{qn}{critical values used to construct joint CI}
#' \item{...}{...}
#'
#' @author Erjia Cui \email{ecui@@umn.edu}, Gabriel Loewinger
#' \email{gloewinger@@gmail.com}
#'
#' @references Cui, E., Leroux, A., Smirnova, E., Crainiceanu, C. (2022). Fast
#' Univariate Inference for Longitudinal Functional Models. \emph{Journal of
#' Computational and Graphical Statistics}, 31(1), 219-230.
#'
#' @export
#'
#' @import lme4
#' @import parallel
#' @import magrittr
#' @import dplyr
#' @import stringr
#' @import mgcv
#' @import refund
#' @importFrom MASS ginv
#' @import cAIC4
#' @importFrom lsei pnnls lsei
#' @import Matrix
#' @importFrom mvtnorm rmvnorm
#' @importFrom Rfast rowMaxs
#' @import progress
#'
#' @examples
#' library(refund)
#'
#' ## random intercept only
#' set.seed(1)
#' DTI_use <- DTI[DTI$ID %in% sample(DTI$ID, 10),]
#' fit_dti <- fui(cca ~ case + visit + sex + (1 | ID),
#'                  data = DTI_use)

fui_conc <- function(
  formula,
  data,
  family = "gaussian",
  var = TRUE,
  analytic = TRUE,
  parallel = FALSE,
  silent = FALSE,
  argvals = NULL,
  nknots_min = NULL,
  nknots_min_cov = 35,
  smooth_method = "GCV.Cp",
  splines = "tp",
  design_mat = FALSE,
  residuals = FALSE,
  G_return = FALSE,
  num_boots = 500,
  boot_type = NULL,
  seed = 1,
  subj_ID = NULL,
  num_cores = 1,
  caic = FALSE,
  REs = FALSE,
  non_neg = 0,
  MoM = 2){

  # If doing parallel computing, set up the number of cores
  if(parallel & !is.integer(num_cores) ) num_cores <- parallel::detectCores() - 1

  # For non-Gaussian family, only do bootstrap inference
  if(family != "gaussian") analytic <- FALSE

  # Organize the input from the model formula
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)

  # stop function if there are column names with "." to avoid issues with covariance G() and H() calculations below
  dep_str <- deparse(model_formula[3])
  if(grepl(".", dep_str, fixed = TRUE)){
    # make sure it isn't just a call to all covariates with "Y ~. "
    dep_str_rm <- substr(dep_str, 3, nchar(dep_str)) # remove first character of parsed formula string and check
    if(grepl(".", dep_str_rm, fixed = TRUE)){
      stop('Remove the character "." from all non-functional covariate names and rerun fui() function
           -For example, change the name "X.1" to "X_1"
           -The string "." *should* be kept in the functional outcome names (e.g., "Y.1" *is* proper naming).')
    }
  }
  rm(dep_str)

  ##############################################################################
  ## Step 1
  ##############################################################################
  if(silent == FALSE) print("Step 1: Fit Massively Univariate Mixed Models")

  # Obtain the dimension of the functional domain
  # indices that start with the outcome name
  out_index <- grep(paste0("^", model_formula[2]), names(data))

  if (length(out_index) != 1) { # functional observations stored in multiple columns
    L <- length(out_index)
  } else { # observations stored as a matrix in one column using the I() function
    L <- ncol(data[,out_index])
  }
  if(analytic & !is.null(argvals) & var)  message("'argvals' argument is currently only supported for bootstrap. `argvals' ignored: fitting model to ALL points on functional domain")
  if(is.null(argvals) | analytic){
    # Set up the functional domain when not specified
    argvals <- 1:L
  }else{
    # when only using a subset of the functional domain
    if(max(argvals) > L) stop("Maximum index specified in argvals is greater than the total number of columns for the functional outcome")
    L <- length(argvals)
  }
  if(family == "gaussian" & analytic & L > 400 & var)   message("Your functional data is dense! Consider subsampling along the functional domain (i.e., reduce columns of outcome matrix) or using bootstrap inference.")

  # Create a matrix to store AICs
  AIC_mat <- matrix(NA, nrow = L, ncol = 2)

  # Fit massively univariate mixed models
  if(parallel == TRUE){
    massmm <- mclapply(
      argvals,
      unimm,
      data = data,
      model_formula = model_formula,
      family = family,
      residuals = residuals,
      caic = caic,
      REs = REs,
      analytic = analytic,
      mc.cores = num_cores
    )
  } else {
    massmm <- mclapply(
      argvals,
      unimm,
      data = data,
      model_formula = model_formula,
      family = family,
      residuals = residuals,
      caic = caic,
      REs = REs,
      analytic = analytic
    )
  }

  # Obtain betaTilde, fixed effects estimates
  betaTilde <- t(do.call(rbind, lapply(massmm, '[[', 1)))
  colnames(betaTilde) <- argvals

  # Obtain residuals, AIC, BIC, and random effects estimates (analytic)
  ## AIC/BIC
  mod_aic <- do.call(c, lapply(massmm, '[[', 3))
  mod_bic <- do.call(c, lapply(massmm, '[[', 4))
  mod_caic <- do.call(c, lapply(massmm, '[[', 6))
  AIC_mat <- cbind(mod_aic, mod_bic, mod_caic)
  colnames(AIC_mat) <- c("AIC", "BIC", "cAIC")
  ## residuals
  resids <- NA
  if(residuals) resids <- suppressMessages(lapply(massmm, '[[', 5) %>% dplyr::bind_cols())
  ## random effects
  if(analytic == TRUE){
    if(REs){
      randEff <- suppressMessages(
        simplify2array(lapply(lapply(massmm, '[[', 7), function(x) as.matrix(x[[1]])))
      )  # will need to change [[1]] to random effect index if multiple REs
    }else{
      randEff <- NULL
    }
    se_mat <- suppressMessages(do.call(cbind, lapply(massmm, '[[', 9) ))
  }else{
    randEff <- se_mat <- NULL
  }

  # Obtain variance estimates of random effects (analytic)
  if(analytic == TRUE){
    var_random <- t(do.call(rbind, lapply(massmm, '[[', 8)))
    sigmaesqHat <- var_random["var.Residual",,drop = FALSE]
    sigmausqHat <- var_random[which(rownames(var_random) != "var.Residual"),,drop = FALSE]
    ## fit a fake model to obtain design matrix
    data$Yl <- unclass(data[,out_index][,1])
    if(family == "gaussian"){
      fit_uni <- suppressMessages(lmer(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                       data = data,
                                       control = lmerControl(optimizer = "bobyqa")))
    }else{
      fit_uni <- suppressMessages(glmer(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                        data = data,
                                        family = family,
                                        control = glmerControl(optimizer = "bobyqa")))
    }
    designmat <- stats::model.matrix(fit_uni) ## model design matrix
    name_random <- as.data.frame(VarCorr(fit_uni))[which(!is.na(as.data.frame(VarCorr(fit_uni))[,3])),3]
    ## names of random effects
    RE_table <- as.data.frame(VarCorr(fit_uni))
    ranEf_grp <- RE_table[,1]
    RE_table <- RE_table[RE_table$grp != "Residual", 1:3]
    ranEf_grp <- ranEf_grp[ranEf_grp != "Residual"]
    ztlist <- sapply(getME(fit_uni, "Ztlist"), function(x) t(x) )
    ## Check if group contains ":" which indicates hierarchical structure and group needs to be specified
    group <- massmm[[1]]$group ## group name in the data
    if(grepl(":", group, fixed = TRUE)){
      if(is.null(subj_ID)){
        group <- str_remove(group, ".*:") # assumes the ID name is to the right of the ":"
      }else{
        group <- subj_ID # use user specified if it exists
      }
    }
    if(is.null(subj_ID)) subj_ID <- group
    randint_flag <- I(length(fit_uni@cnms) == 1 & length(fit_uni@cnms[[group]]) == 1 & fit_uni@cnms[[group]][1] == "(Intercept)")
    rm(fit_uni)
  }


  ##############################################################################
  ## Step 2
  ##############################################################################
  if(silent == FALSE) print("Step 2: Smoothing")

  # Penalized splines smoothing and extract components (analytic)
  nknots <- min(round(L/2), nknots_min) ## number of knots for regression coefficients
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov) ## number of knots for covariance matrix
  nknots_fpca <- min(round(L/2), 35)
  ## smoothing parameter, spline basis, penalty matrix (analytic)
  if(analytic == TRUE){
    p <- nrow(betaTilde) ## number of fixed effects parameters
    betaHat <- matrix(NA, nrow = p, ncol = L)
    lambda <- rep(NA, p)
    for(r in 1:p){
      fit_smooth <- gam(betaTilde[r,] ~ s(argvals, bs = splines, k = (nknots + 1)), method = smooth_method)
      betaHat[r,] <- fit_smooth$fitted.values
      lambda[r] <- fit_smooth$sp ## smoothing parameter
    }
    sm <- smoothCon(s(argvals, bs = splines, k = (nknots + 1)),
                    data=data.frame(argvals=argvals),
                    absorb.cons=TRUE)
    S <- sm[[1]]$S[[1]] ## penalty matrix
    B <- sm[[1]]$X ## basis functions
    rm(fit_smooth, sm)
  }else{
    betaHat <- t(apply(betaTilde, 1, function(x) gam(x ~ s(argvals, bs = splines, k = (nknots + 1)), method = smooth_method)$fitted.values))
  }
  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- 1:L


  ##############################################################################
  ## Step 3
  ##############################################################################
  if(var == TRUE){ ## skip the step when var = FALSE

    if(analytic == TRUE){

      ##########################################################################
      ## Analytic Inference
      ##########################################################################
      if(silent == FALSE) print("Step 3: Inference (Analytic)")

      ##########################################################################
      ## Step 3.1
      ##########################################################################
      if(silent == FALSE) print("Step 3.1: Preparation")

      # Fill in missing values of the raw data using FPCA

      if(length(which(is.na(data[,out_index]))) != 0){
        message(paste("Imputing", length(out_index), "values in functional response with longitudinal functional PCA" ))
        if(length(out_index) != 1){
          tmp <- as.matrix(data[,out_index])
          tmp[which(is.na(tmp))] <-
            suppressWarnings( refund::fpca.face(tmp,
                                                argvals = argvals,
                                                knots = nknots_fpca)$Yhat[which(is.na(tmp))] )
          data[,out_index] <- tmp
        }else{
          data[,out_index][which(is.na(data[,out_index]))] <-
            suppressWarnings( refund::fpca.face(as.matrix(data[,out_index]),
                                                argvals = argvals,
                                                knots = nknots_fpca)$Yhat[which(is.na(data[,out_index]))] )
        }
      }

      # Derive variance estimates of random components: H(s), R(s)
      HHat <- t(apply(sigmausqHat, 1, function(b) stats::smooth.spline(x = argvals, y = b)$y))
      ind_var <- which(grepl("var", rownames(HHat)) == TRUE) ## index of variance
      HHat[ind_var,][which(HHat[ind_var,] < 0)] <- 0
      RHat <- t(apply(sigmaesqHat, 1, function(b) stats::smooth.spline(x = argvals, y = b)$y))
      RHat[which(RHat < 0)] <- 0

      ### Create a function that estimates covariance G for general cases
      #### Create a subfunction that estimates non-negative terms on diagonal
      cov.nnls <- function(data_cov, RE_table, idx_lst, d_temp, designmat,
                           betaHat, GTilde, non_neg = 0, silent = TRUE){

        if(non_neg == 1){
          if(silent == FALSE) print("Step 3.1.2: NNLS 1")

          # put constraints on EVERY coef corresponding to columns for one random effect
          ncol_Z <- ncol(data_cov$Z)
          var_term_idx <- which(is.na(RE_table$var2)) # find variance terms (that need non-negativity)

          idx_start_end <- rep(NA, ncol_Z ) # initial indices
          non_negIdx <- do.call(c, lapply(var_term_idx, function(xx) idx_lst[[xx]]   ) ) # concatenate all the terms that correspond to non-negative constraints
          kk <- ncol_Z - length(non_negIdx) # indices of unconstrained coefficients
          idx_start_end[1:kk] <- seq(1,ncol_Z)[-non_negIdx] # put remaining terms after
          idx_start_end[(kk+1):ncol_Z] <- non_negIdx # put unconstrained terms first
          XX <- as.matrix(data_cov$Z[,idx_start_end]) # reorder to put covariates associated with non-negative coefs in first columns
          bHat <- rep(NA, dim(GTilde)[1])

          for(i in 1:L){
            YYi <- (d_temp[,i] - designmat %*% betaHat[,i])^2 # outcome of product of residuals
            bHat[idx_start_end] <- lsei::pnnls(a = XX, b = YYi, k = kk)$x # nnls
            GTilde[,i,i]  <- sapply(idx_lst,  function(x) mean(bHat[x]) ) # average over coefficients corresponding to same random effect term
          }

        }else if(non_neg == 2){
          if(silent == FALSE) print("Step 3.1.2: NNLS 2")

          # put constraints on AVERAGE over coefs corresponding to columns for one random effect
          ncol_Z <- ncol(data_cov$Z)
          ff <- rep(0, nrow(RE_table)) # non-negativity vector
          eMat <- matrix(0, nrow = nrow(RE_table), ncol = ncol_Z) # constraint vector to be matrix below ( initially make all 0s so we do not place constraints on terms that can be negative)
          var_term_idx <- which(is.na(RE_table$var2)) # find variance terms (that need non-negativity)

          for(ii in var_term_idx){
            eMat[ii, idx_lst[[ii]] ] <- 1 # use these sum (and/or average) to enforce constraint on average
          }

          XX <- as.matrix(data_cov$Z) # reorder to put covariates associated with non-negative coefs in first columns and then transpose design mat for package

          for(i in 1:L){
            YYi <- (d_temp[,i] - designmat %*% betaHat[,i])^2 # outcome of product of residuals
            bHat <- lsei::lsei(a = XX, b = YYi, e = eMat, f = ff) # nnls -- allows tiny negative values due to error
            GTilde[,i,i]  <- sapply(idx_lst,  function(x) mean(bHat[x]) ) # average over coefficients corresponding to same random effect term
          }

        }

        return(GTilde)
      }

      ### Create a function that takes in design matrix and names to calculate
      ### and produce design matrix and indices for summing for G(s_1, s_2)
      ### covariance matrix OLS regression
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
              Z[[j]] <- Z_lst[[ zlst_idx ]]
            }else{
              # not an intercept term (corresponds to random slope so requires squaring elements) -- see below for why we can square instead of doing actual element-wise produt (all other product terms are zero-ed out)
              Z[[j]] <- (Z_lst[[ zlst_idx ]])^2
            }

          }else{
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

          if(ID_flag){
            Z[[j]] <- matrix( rowSums(Z[[j]]), ncol = 1)
            colnames(Z[[j]]) <- paste0(RE_table$grp[j], "_", RE_table$var1[j], v2) # name column
            idx_vec[j] <- 1
          }else{
            idx_vec[j] <- ncol(Z[[j]])         # number of columns in matrix
            colnames(Z[[j]]) <- paste0(RE_table$grp[j], "_", RE_table$var1[j], v2, "_", 1:ncol(Z[[j]]))   # name columns
          }

        }

        idx_vec <- c(0, cumsum(idx_vec))
        idx_lst <- sapply( seq_along(1:nrow(RE_table)),  function(x)   (idx_vec[x] + 1):(idx_vec[x+1]) ) # column indices
        Z <- do.call(cbind, Z) # concatenate

        return( list(Z = Z,
                     Z_orig = Z_orig,
                     idx_lst = idx_lst,
                     idx_vec = idx_vec[-1])
        )

      }

      # TODO: Place this in its own function
      # Also needs to call `select_knots` and `pspline_setting` w/in the pkg

      message("Testing pspline_setting and select_knots")
      # altered fbps() function from refund to increase GCV speed (b/c lambda1=lambda2 for cov matrices b/c they're symmetric)
      fbps_cov <- function(data, subj=NULL,covariates = NULL, knots=35, knots.option="equally-spaced",
                       periodicity = c(FALSE,FALSE), p=3,m=2,lambda=NULL,
                       selection = "GCV",
                       search.grid = T, search.length = 100, method="L-BFGS-B",
                       lower= -20, upper=20, control=NULL){

        # return a smoothed matrix using fbps

        # data: a matrix
        # covariates: the list of data points for each dimension
        # knots: to specify either the number/numbers of  knots  or the vector/vectors of knots for each dimension; defaults to 35
        # p: the degrees of B-splines
        # m: the order of difference penalty
        # lambda: the user-selected smoothing parameters
        # lscv: for leave-one-subject-out cross validation, the columns are subjects
        # method: see "optim"
        # lower, upper, control: see "optim"

        ## data dimension
        data_dim = dim(data)
        n1 = data_dim[1]
        n2 = data_dim[2]

        ## subject ID
        if(is.null(subj)) subj = 1:n2
        subj_unique = unique(subj)
        I = length(subj_unique)
        ## covariates for the two axis
        if(!is.list(covariates)) {

          x=(1:n1)/n1-1/2/n1; ## if NULL, assume equally distributed
          z = (1:n2)/n2-1/2/n2
        }
        if(is.list(covariates)){

          x = covariates[[1]]
          z = covariates[[2]]
        }

        ## B-spline basis setting
        p1 = rep(p,2)[1]
        p2 = rep(p,2)[2]
        m1 = rep(m,2)[1]
        m2 = rep(m,2)[2]

        ## knots
        if(!is.list(knots)){
          K1 = rep(knots,2)[1]
          xknots = select_knots(x,knots=K1,option=knots.option)
          K2 = rep(knots,2)[2]
          zknots = select_knots(z,knots=K2,option=knots.option)
        }

        if(is.list(knots)){

          xknots = knots[[1]]
          K1 = length(xknots)-1
          knots_left <- 2*xknots[1]-xknots[p1:1+1]
          knots_right <- 2*xknots[K1] - xknots[K1-(1:p1)]
          xknots <- c(knots_left,xknots,knots_right)

          zknots= knots[[2]]
          K2 = length(zknots)-1
          knots_left <- 2*zknots[1]- zknots[p2:1+1]
          knots_right <- 2*zknots[K2] - zknots[K2-(1:p2)]
          zknots <- c(knots_left,zknots,knots_right)
        }
        #######################################################################
        Y = data

        ###################   precalculation for fbps smoothing  ##############

        List = pspline_setting(x,xknots,p1,m1,periodicity[1])
        A1 = List$A
        B1 = List$B
        Bt1 = Matrix(t(as.matrix(B1)))
        s1 = List$s
        Sigi1_sqrt = List$Sigi_sqrt
        U1 = List$U
        A01 = Sigi1_sqrt%*%U1
        c1 = length(s1)

        List = pspline.setting(z,zknots,p2,m2,periodicity[2])
        A2 = List$A
        B2 = List$B
        Bt2 = Matrix(t(as.matrix(B2)))
        s2 = List$s
        Sigi2_sqrt = List$Sigi.sqrt
        U2 = List$U
        A02 = Sigi2_sqrt%*%U2
        c2 = length(s2)
        #################select optimal penalty ################################

        tr <-function(A){ return(sum(diag(A)))} ## the trace of a square matrix

        Ytilde = Bt1%*%(Y%*%B2)
        Ytilde = t(A01)%*%Ytilde%*%A02
        Y_sum = sum(Y^2)
        ytilde = as.vector(Ytilde)
        if(selection=="iGCV"){

          KH = function(A,B){
            C = matrix(0,dim(A)[1],dim(A)[2]*dim(B)[2])
            for(i in 1:dim(A)[1])
              C[i,] = Matrix::kronecker(A[i,],B[i,])
            return(C)
          }
          G = rep(0,I)
          Ybar = matrix(0,c1,I)
          C = matrix(0,c2,I)

          Y2 = Bt1%*%Y
          Y2 = matrix(t(A01)%*%Y2,c1,n2)
          for(i in 1:I){
            sel = (1:n2)[subj==subj_unique[i]]
            len = length(sel)
            G[i] = len
            Ybar[,i] = as.vector(matrix(Y2[,sel],ncol=len)%*%rep(1,len))
            C[,i] = as.vector(t(matrix(A2[sel,],nrow=len))%*%rep(1,len))
          }
          g1 = diag(Ybar%*%diag(G)%*%t(Ybar))
          g2 = diag(Ybar%*%t(Ybar))
          g3 = ytilde*as.vector(Ybar%*%diag(G)%*%t(C))
          g4 = ytilde*as.vector(Ybar%*%t(C))
          g5 = diag(C%*%diag(G)%*%t(C))
          g6 = diag(C%*%t(C))
          #cat("Processing completed\n")
        }

        fbps_gcv =function(x){

          lambda=exp(x)
          ## two lambda's are the same
          if(length(lambda)==1)
          {
            lambda1 = lambda
            lambda2 = lambda
          }
          ## two lambda's are different
          if(length(lambda)==2){
            lambda1=lambda[1]
            lambda2=lambda[2]
          }

          sigma2 = 1/(1+lambda2*s2)
          sigma1 = 1/(1+lambda1*s1)
          sigma2_sum = sum(sigma2)
          sigma1_sum = sum(sigma1)
          sigma = Matrix::kronecker(sigma2,sigma1)
          sigma.2 = Matrix::kronecker(sqrt(sigma2),sqrt(sigma1))

          if(selection=="iGCV"){
            d = 1/(1-(1+lambda1*s1)/sigma2_sum*n2)
            gcv = sum((ytilde*sigma)^2) - 2*sum((ytilde*sigma.2)^2)
            gcv = gcv + sum(d^2*g1) - 2*sum(d*g2)
            gcv = gcv - 2*sum(g3*Matrix::kronecker(sigma2,sigma1*d^2))
            gcv = gcv + 4*sum(g4*Matrix::kronecker(sigma2,sigma1*d))
            gcv = gcv + sum(ytilde^2*Matrix::kronecker(sigma2^2*g5,sigma1^2*d^2))
            gcv = gcv - 2*sum(ytilde^2*Matrix::kronecker(sigma2^2*g6,sigma1^2*d))

          }

          if(selection=="GCV") {

            gcv = sum((ytilde*sigma)^2) - 2*sum((ytilde*sigma.2)^2)
            gcv = Y_sum + gcv
            trc = sigma2_sum*sigma1_sum
            gcv = gcv/(1-trc/(n1*n2))^2
          }
          return(gcv)
        }

        fbps_est =function(x){

          lambda=exp(x)
          ## two lambda's are the same
          if(length(lambda)==1)
          {
            lambda1 = lambda
            lambda2 = lambda
          }
          ## two lambda's are different
          if(length(lambda)==2){
            lambda1=lambda[1]
            lambda2=lambda[2]
          }

          sigma2 = 1/(1+lambda2*s2)
          sigma1 = 1/(1+lambda1*s1)
          sigma2_sum = sum(sigma2)
          sigma1_sum = sum(sigma1)
          sigma = Matrix::kronecker(sigma2,sigma1)
          sigma.2 = Matrix::kronecker(sqrt(sigma2),sqrt(sigma1))

          Theta = A01%*%diag(sigma1)%*%Ytilde
          Theta = as.matrix(Theta%*%diag(sigma2)%*%t(A02))
          Yhat = as.matrix(as.matrix(B1%*%Theta)%*%Bt2)
          result=list(lambda=c(lambda1,lambda2),Yhat=Yhat,Theta=Theta,
                      setting = list(x = list(knots = xknots, p = p1, m = m1),
                                     z = list(knots = zknots, p = p2, m = m2)))
          class(result) ="fbps"
          return(result)
        }

        if(is.null(lambda)){

          if(search.grid ==T){
            lower2 <- lower1 <- lower[1]
            upper2 <- upper1 <- upper[1]
            search.length2 <- search.length1 <- search.length[1]
            if(length(lower)==2) lower2 <- lower[2]
            if(length(upper)==2) upper2 <- upper[2]
            if(length(search.length)==2) search.length2 <- search.length[2]

            Lambda1 = seq(lower1,upper1,length = search.length1)
            lambda.length1 = length(Lambda1)

            GCV = rep(0, lambda.length1)
            for(j in 1:lambda.length1){
              GCV[j] = fbps_gcv(c(Lambda1[j], Lambda1[j]))
            }

            location = which.min(GCV)[1]
            j0 = location%%lambda.length1
            if(j0==0) j0 = lambda.length1
            k0 = (location-j0)/lambda.length1+1
            lambda = exp(c(Lambda1[j0], Lambda1[j0]))
          } ## end of search.grid

          if(search.grid == F){
            fit = stats::optim(0,fbps_gcv,method=method,control=control,
                        lower=lower[1],upper=upper[1])

            fit = stats::optim(c(fit$par,fit$par),fbps_gcv,method=method,control=control,
                        lower=lower[1],upper=upper[1])
            if(fit$convergence>0) {
              expression = paste("Smoothing failed! The code is:",fit$convergence)
              print(expression)
            }
            lambda = exp(fit$par)
          } ## end of optim

        } ## end of finding smoothing parameters
        lambda = rep(lambda,2)[1:2]

        return(fbps_est(log(lambda)))
      }


      ## Calculate Method of Moments estimator for G() with potential NNLS correction for diagonals (for variance terms)
      if(randint_flag) {

        message("Testing integration of G_estimate_randint")

        GTilde <- G_estimate_randint(
          data = data, out_index = out_index,
          designmat = designmat, betaHat = betaHat,
          silent = silent
        )

        message("G_estimate_randint successsful.")

        if(silent == FALSE)
          print("Step 3.1.2: Smooth G")

        diag(GTilde) <- HHat[1,] # L x L matrix
        GHat <- fbps_cov(GTilde, search.length = 100, knots = nknots_cov)$Yhat ## fast bivariate smoother # nknots_min
        diag(GHat)[which(diag(GHat) < 0)] <- diag(GTilde)[which(diag(GHat) < 0)]

      } else {
        ## if more random effects than random intercept only
        if(silent == FALSE)
          print("Step 3.1.1: Preparation B")

        # generate design matrix for G(s_1, s_2)
        # Method of Moments Linear Regression calculation

        data_cov <- G_generate(data = data,
                               Z_lst = ztlist,
                               RE_table = RE_table,
                               ID = subj_ID)

        ## Method of Moments estimator for G() with potential NNLS correction
        # for diagonals (for variance terms)

        GTilde <- G_estimate(data = data,
                             out_index = out_index,
                             data_cov = data_cov,
                             ztlist = ztlist,
                             designmat = designmat,
                             betaHat = betaHat,
                             HHat = HHat,
                             RE_table = RE_table,
                             non_neg = non_neg,
                             MoM = MoM,
                             silent = silent) # method of moments estimator

        if(silent == FALSE) print("Step 3.1.2: Smooth G")

        ## smooth GHat
        GHat <- GTilde
        for(r in 1:nrow(HHat)){
          diag(GTilde[r,,]) <- HHat[r,]
          GHat[r,,] <- fbps_cov(GTilde[r,,], search.length = 100, knots = nknots_cov)$Yhat # , knots = nknots_min
        }
        diag(GHat[r,,])[which(diag(GHat[r,,]) < 0)] <- diag(GTilde[r,,])[which(diag(GHat[r,,]) < 0)]
      }

      # For return values below
      if(!design_mat ){
        ztlist <- NULL
        idx_vec_zlst <- NULL
      }

      ##########################################################################
      ## Step 3.2
      ##########################################################################
      if(silent == FALSE) print("Step 3.2: First step")

      # Calculate the intra-location variance of betaTilde: Var(betaTilde(s))

      # fast block diagonal generator taken from Matrix package examples
      bdiag_m <- function(lmat) {
        ## Copyright (C) 2016 Martin Maechler, ETH Zurich
        if(!length(lmat)) return(methods::new("dgCMatrix"))
        stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
                  (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
                  all(vapply(lmat, dim, integer(2)) == k)) # all of them
        N <- length(lmat)
        if(N * k > .Machine$integer.max)
          stop("resulting matrix too large; would be  M x M, with M=", N*k)
        M <- as.integer(N * k)
        ## result: an   M x M  matrix
        methods::new("dgCMatrix", Dim = c(M,M),
            ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
            i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
            p = k * 0L:M,
            x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
      }

      ### Create a function that organizes covariance matrices correctly based on RE table
      # use this to find indices to feed into  cov_organize  function above
      cov_organize_start <- function(cov_vec){
        # assumes each set of cov for 2 preceeding variance terms ofrandom effects
        # MAY NEED TO BE UPDATED FOR MORE COMPLICATED RANDOM EFFECT STRUCTURES
        # the order it spits out is the order given by cov_vec (we use (sort below))

        nm <- names(cov_vec) # names of the terms to determine which are variances and which are covariances
        cov_idx <- which( grepl("cov", nm, fixed = TRUE) ) # find terms that are covariances
        covs <- length(cov_idx) # number of covariances

        # organize based on number of covariance terms
        if(covs == 0){
          # if no covariance terms (just simple diagonal matrix)
          var_nm <- groupings <- groups_t <- g_idx_list <- v_list <- NULL
          v_list_template <- diag(1:length(cov_vec))
          v_list_template <- apply(v_list_template, 1, function(x) ifelse(x == 0, max(v_list_template) + 1, x) ) # replace 0s with corresponding index

        }else{
          # mix of covariance terms and non-covariance terms
          # variance terms
          var_nm <- nm[sapply(nm, function(x) unlist(strsplit(x, split='.', fixed=TRUE))[1] == "var"  )]
          groupings <- unique( sapply(var_nm, function(x) unlist(strsplit(x, split='.', fixed=TRUE))[2] ) )
          g_idx_list <- mat_lst <- vector(length = length(groupings), "list")
          cnt <- 0
          groups_t <- sapply(var_nm, function(x) unlist(strsplit(x, split='.', fixed=TRUE))[2] ) # grouping variable for each name
          v_list <- vector(length = length(groupings), "list")

          # iterate through groupings and make sub-matrices
          for(trm in groupings){
            cnt <- cnt + 1
            # find current grouping (e.g., "id" or "id:session")
            g_idx <- g_idx_list[[cnt]] <- names(which( groups_t == trm  )) # current grouping
            # if this is not just a variance term (i.e., cov terms too)
            if(length(g_idx) > 1){
              m <- diag(cov_vec[g_idx]) # make diagonal matrix with variance terms
              v_list[[cnt]] <- matrix(NA, ncol = ncol(m), nrow = ncol(m))

              # name after individual random effect variables (covariates)
              trm_var_name <- sapply(g_idx, function(x) unlist(strsplit(x, split='.', fixed=TRUE))[3] )
              nm1 <- paste0("cov.", trm, ".")
              for(ll in 1:ncol(m)){
                for(jj in seq(1,ncol(m) )[-ll] ){
                  # names that are covariance between relevant variables
                  v1 <- which(nm == paste0(nm1, trm_var_name[ll], ".", trm_var_name[jj]) )
                  v2 <- which(nm == paste0(nm1, trm_var_name[jj], ".", trm_var_name[ll]) )
                  v_list[[cnt]][ll,jj] <- c(v1, v2)
                }
              }
              v_list[[cnt]][is.na(v_list[[cnt]])] <- sapply(g_idx, function(x) which(nm == x) )
            }
          }

          v_list_template <- bdiag_m(v_list) # diagonal matrix of indices for template
          v_list_template <- apply(v_list_template, 1, function(x) ifelse(x == 0, max(v_list_template) + 1, x) ) # replace 0s with corresponding index (just highest + 1)
        }

        return( list(nm=nm,
                     cov_idx=cov_idx,
                     covs=covs,
                     groupings=groupings,
                     g_idx_list=g_idx_list,
                     v_list = v_list,
                     v_list_template = v_list_template) )
      }

      ### Create a function that trims eigenvalues and return PSD matrix, V is a n x n matrix
      eigenval_trim <- function(V){
        edcomp <- eigen(V, symmetric = TRUE) ## trim non-positive eigenvalues to ensure positive semidefinite
        eigen.positive <- which(edcomp$values > 0)
        q=ncol(V)

        if(length(eigen.positive) == q){
          # nothing needed here because matrix is already PSD
          return(V)
        }else if(length(eigen.positive) == 0){
          return(tcrossprod(edcomp$vectors[,1]) * edcomp$values[1])
        }else if(length(eigen.positive) == 1){
          return(tcrossprod(as.vector(edcomp$vectors[,1])) * as.numeric(edcomp$values[1]))
        }else{
          # sum of outerproducts of eigenvectors scaled by eigenvalues for all positive eigenvalues
          return(matrix(edcomp$vectors[,eigen.positive] %*% tcrossprod(Diagonal(x=edcomp$values[eigen.positive]), edcomp$vectors[,eigen.positive]), ncol = q))
        }
      }

      ## Obtain the corresponding rows of each subject
      obs.ind <- list()
      ID.number <- unique(data[,group])
      for(id in ID.number){
        obs.ind[[as.character(id)]] <- which(data[,group] == id)
      }
      ## Concatenate vector of 1s to Z because used that way below
      HHat_trim <- NA
      if(!randint_flag){
        Z <- data_cov$Z_orig
        qq <- ncol(Z)
        HHat_trim <- array(NA, c(qq, qq, L)) # array for Hhat
      }
      ## Create var.beta.tilde.theo to store variance of betaTilde
      var.beta.tilde.theo <- array(NA, dim = c(p,p,L))
      ## Create XTVinvZ_all to store all XTVinvZ used in the covariance calculation
      XTVinvZ_all <- vector(length = L, "list")
      ## arbitrarily start find indices
      resStart <- cov_organize_start(HHat[,1])
      res_template <- resStart$v_list_template # index template
      template_cols <- ncol(res_template)
      ## Calculate Var(betaTilde) for each location
      for(s in 1:L){
        V.subj.inv <- c()
        ## we first do inverse of each block matrix then combine them
        if(!randint_flag){
          cov.trimmed <- eigenval_trim( matrix(c(HHat[,s], 0)[res_template], template_cols) )
          HHat_trim[,,s] <- cov.trimmed
        }

        XTVinvX <- matrix(0, nrow = p, ncol = p) # store XT * Vinv * X
        XTVinvZ_i <- vector(length = length(ID.number), "list") # store XT * Vinv * Z
        for(id in ID.number){ ## iterate for each subject
          subj.ind <- obs.ind[[as.character(id)]]
          if(randint_flag){
            Ji <- length(subj.ind)
            V.subj <- matrix(HHat[1,s], nrow = Ji, ncol = Ji) + diag(RHat[s], Ji)
          }else{
            V.subj <- Z[subj.ind,,drop = FALSE] %*% tcrossprod( cov.trimmed, Z[subj.ind,,drop = FALSE] ) +
              diag(RHat[s], length(subj.ind))
          }

          V.subj.inv <- solve(V.subj)

          XTVinvX <- XTVinvX + crossprod(matrix(designmat[subj.ind,], ncol = p),
                                         V.subj.inv) %*% matrix(designmat[subj.ind,], ncol = p)
          if(randint_flag){
            XTVinvZ_i[[as.character(id)]] <- crossprod(matrix(designmat[subj.ind,], ncol = p),
                                                       V.subj.inv) %*% matrix(1, nrow = Ji, ncol = 1)
          }else{
            XTVinvZ_i[[as.character(id)]] <- crossprod(matrix(designmat[subj.ind,], ncol = p),
                                                       V.subj.inv) %*% Z[subj.ind,,drop = FALSE]
          }
        }

        var.beta.tilde.theo[,,s] <- solve(XTVinvX) # variance of betaTilde
        XTVinvZ_all[[s]] <- XTVinvZ_i # store XTVinvZ for each subject at each location

      }
      suppressWarnings(rm(V.subj.inv, cov.trimmed, XTVinvZ_i, resStart, res_template, template_cols))

      # Calculate the inter-location covariance of betaTilde: Cov(betaTilde(s_1), betaTilde(s_2))
      ## Create cov.beta.tilde.theo to store covariance of betaTilde
      cov.beta.tilde.theo <- array(NA, dim = c(p,p,L,L))
      if(randint_flag){
        resStart <- cov_organize_start(GHat[1,2]) # arbitrarily start
      }else{
        resStart <- cov_organize_start(GHat[,1,2]) # arbitrarily start
      }

      res_template <- resStart$v_list_template # index template
      template_cols <- ncol(res_template)
      ## Calculate Cov(betaTilde) for each pair of location
      for(i in 1:L){
        for(j in i:L){
          V.cov.subj <- list()
          tmp <- matrix(0, nrow = p, ncol = p) ## store intermediate part
          if(randint_flag){
            G_use <- GHat[i,j]
          }else{
            G_use <- eigenval_trim( matrix(c(GHat[,i,j], 0)[res_template], template_cols) )
          }

          for(id in ID.number){
            tmp <- tmp + XTVinvZ_all[[i]][[as.character(id)]] %*% tcrossprod(G_use, XTVinvZ_all[[j]][[as.character(id)]])
          }

          ## Calculate covariance using XTVinvX and tmp to save memory
          cov.beta.tilde.theo[,,i,j] <- var.beta.tilde.theo[,,i] %*% tmp %*% var.beta.tilde.theo[,,j]

        }
      }
      suppressWarnings(rm(V.subj, V.cov.subj, Z, XTVinvZ_all, resStart, res_template, template_cols))


      ##########################################################################
      ## Step 3.3
      ##########################################################################
      if(silent == FALSE) print("Step 3.3: Second step")

      # Do an intermediate step for covariance estimate
      var.beta.tilde.s <- array(NA, dim = c(L,L,p))
      for(j in 1:p){
        for(r in 1:L){
          for(t in 1:L){
            if(t == r){
              var.beta.tilde.s[r,t,j] <- var.beta.tilde.theo[j,j,r]
            }else{
              var.beta.tilde.s[r,t,j] <- cov.beta.tilde.theo[j,j,min(r,t),max(r,t)]
            }
          }
        }
      }

      # Calculate the inter-location covariance of betaHat: Cov(betaHat(s_1), betaHat(s_2))
      var.beta.hat <- array(NA, dim = c(L,L,p))
      for(r in 1:p){
        M <- B %*% tcrossprod(solve(crossprod(B) + lambda[r]*S), B) + matrix(1/L, nrow = L, ncol = L)
        var.raw <- M %*% tcrossprod(var.beta.tilde.s[,,r], M)
        var.beta.hat[,,r] <- eigenval_trim(var.raw) ## trim eigenvalues to make final variance matrix PSD
      }
      betaHat.var <- var.beta.hat ## final estimate

      # Obtain qn to construct joint CI
      qn <- rep(0, length = nrow(betaHat))
      N <- 10000 ## sample size in simulation-based approach
      zero_vec <- rep(0, length(betaHat[1,]))
      set.seed(seed)

      for(i in 1:length(qn)){
        Sigma <- betaHat.var[,,i]
        sqrt_Sigma <- sqrt(diag(Sigma))
        S_scl <- Matrix::Diagonal(x = 1/sqrt_Sigma)
        Sigma <- as.matrix(S_scl %*% Sigma %*% S_scl)
        x_sample <- abs(mvtnorm::rmvnorm(N, zero_vec, Sigma)) #x_sample <- abs(FastGP::rcpp_rmvnorm_stable(N, Sigma, zero_vec))
        un <- Rfast::rowMaxs(x_sample, value = TRUE)
        qn[i] <- stats::quantile(un, 0.95)
      }

      # Decide whether to return design matrix or just set it to NULL
      if(!design_mat) designmat <- NA
      if(!silent)  message("Complete! \n -Use plot_fui() function to plot estimates. \n -For more information, run the command:  ?plot_fui")

      return(list(betaHat = betaHat, betaHat.var = betaHat.var, qn = qn, aic = AIC_mat,
                  betaTilde = betaTilde, var_random = var_random, designmat = designmat, residuals = resids,
                  H = HHat_trim, R = RHat, G = GTilde, GHat = GHat, Z = ztlist,
                  argvals = argvals,
                  randEff = randEff, se_mat = se_mat))

    } else {

      ##########################################################################
      ## Bootstrap Inference
      ##########################################################################
      if(silent == FALSE) print("Step 3: Inference (Bootstrap)")

      # Check to see if group contains ":" which indicates hierarchical structure and group needs to be specified
      group <- massmm[[1]]$group
      if(grepl(":", group, fixed = TRUE)){
        if(is.null(subj_ID)){
          group <- str_remove(group, ".*:") # assumes the ID name is to the right of the ":"
        }else if(!is.null(subj_ID)){
          group <- subj_ID # use user specified if it exists
        }else{
          message("You must specify the argument: ID")
        }
      }
      ID.number <- unique(data[,group])
      idx_perm <- t( replicate(num_boots, sample.int(length(ID.number), length(ID.number), replace = TRUE)) )
      B <- num_boots
      betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))

      if(is.null(boot_type)){
        # default bootstrap type if not specified
        if(family == "gaussian"){
          boot_type <- ifelse(length(ID.number) <= 10, "reb", "wild")
        }else{
          boot_type <- "cluster"
        }
      }

      if(family != "gaussian" & boot_type %in% c("wild", "reb")){
        stop('Non-gaussian outcomes only supported for some bootstrap procedures. \n Set argument `boot_type` to one of the following: "parametric", "semiparametric", "cluster", "case", "residual"')
      }
      message(paste("Bootstrapping Procedure:", as.character(boot_type)))
      # original way
      if(boot_type == "cluster"){
      # Do bootstrap
        pb <- progress_bar$new(total = B)
        for(boots in 1:B){
          pb$tick()
          sample.ind <- idx_perm[boots,] # take one of the randomly sampled (and unique) combinations
          dat.ind <- new_ids <- vector(length = length(sample.ind), "list")
          for(ii in 1:length(sample.ind)){
            dat.ind[[ii]] <- which(data[,group] == ID.number[sample.ind[ii]])
            new_ids[[ii]] <- rep(ii, length(dat.ind[[ii]]))  # subj_b is now the pseudo_id
          }
          dat.ind <- do.call(c, dat.ind)
          new_ids <- do.call(c, new_ids)
          df2 <- data[dat.ind,] # copy dataset with subset of rows we want
          df2[,subj_ID] <- new_ids # replace old IDs with new IDs

          fit_boot <- fui(formula = formula,
                          data = df2,
                          family = family,
                          argvals = argvals,
                          var = FALSE,
                          parallel = FALSE,
                          silent = TRUE,
                          nknots_min = nknots_min,
                          nknots_min_cov = nknots_min_cov,
                          smooth_method = smooth_method,
                          splines = splines,
                          residuals = residuals,
                          subj_ID = subj_ID,
                          num_cores = num_cores,
                          REs = FALSE)

          betaHat_boot[,,boots] <- fit_boot$betaHat
        }
        rm(fit_boot, df2, dat.ind, new_ids)
      }else{

        # lmeresampler() way
        B <- num_boots # use original amount do not constrain by number of unique resampled types here because we cannot construct rows to resample
        betaHat_boot <- betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
        if(!silent)   print("Step 3.1: Bootstrap resampling-") # , as.character(boot_type)
        pb <- progress_bar$new(total = L)
        for(l in 1:L){
          pb$tick()
          data$Yl <- unclass(data[,out_index][,argvals[l]])
          fit_uni <- suppressMessages(lmer(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                           data = data,
                                           control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5000))))
          # if(l == 1){
          #   fit_uni <- suppressMessages(lmer(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
          #                                  data = data,
          #                                  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5000))))
          #   fit_orig <- fit_uni # use this in refit() so we aren't constantly updating a refit model
          # }else{
          #   fit_uni <- suppressMessages( lme4::refit(fit_orig, newresp = unclass(data[,out_index][,l]),
          #                                            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5000))) )
          # }

          set.seed(seed) # set seed to make sure bootstrap replicate (draws) are correlated across functional domains

          if(boot_type == "residual"){
            boot_sample <- lmeresampler::bootstrap(model = fit_uni,
                                                   B = B,
                                                   type = boot_type,
                                                   rbootnoise = 0.0001)$replicates # for residual bootstrap to avoid singularity problems
            betaTilde_boot[,l,] <- t(as.matrix(boot_sample[,1:nrow(betaHat)]))
          }else if(boot_type %in% c("wild", "reb", "case") ){
            # for case
            flist <- lme4::getME(fit_uni, "flist")
            re_names <- names(flist)
            clusters_vec <- c(rev(re_names), ".id")
            resample_vec <- c(TRUE, rep(FALSE, length(clusters_vec) - 1)) # for case bootstrap, we only resample at first (subject level) because doesn't make sense to resample within-cluster for longitudinal data
            boot_sample <- lmeresampler::bootstrap(model = fit_uni,
                                                   B = B,
                                                   type = boot_type,
                                                   resample = resample_vec, # only matters for type = "case"
                                                   hccme = "hc2", # wild bootstrap
                                                   aux.dist = "mammen", # wild bootstrap
                                                   reb_type = 0)$replicates # for reb bootstrap only
            betaTilde_boot[,l,] <- t(as.matrix(boot_sample[,1:nrow(betaHat)]))
          }else{
            use.u <- ifelse(boot_type == "semiparametric", TRUE, FALSE)
            betaTilde_boot[,l,] <- t( lme4::bootMer(x = fit_uni, FUN = function(.){fixef(.)},
                                                    nsim = B, seed = seed, type = boot_type,
                                                    use.u = use.u)$t )
          }
        }

        suppressWarnings(rm(boot_sample, fit_uni))
        # smooth across functional domain
        if(!silent)   print("Step 3.2: Smooth Bootstrap estimates")
        for(b in 1:B){
          betaHat_boot[,,b] <- t(apply(betaTilde_boot[,,b], 1, function(x) gam(x ~ s(argvals, bs = splines, k = (nknots + 1)), method = smooth_method)$fitted.values))
        }

        rm(betaTilde_boot)
      }

      # Obtain bootstrap variance
      betaHat.var <- array(NA, dim = c(L,L,nrow(betaHat)))
      for(r in 1:nrow(betaHat)){
        betaHat.var[,,r] <- 1.2*var(t(betaHat_boot[r,,])) ## account for within-subject correlation
      }

      # Obtain qn to construct joint CI using the fast approach
      if(!silent)   print("Step 3.3")
      qn <- rep(0, length = nrow(betaHat))
      N <- 10000 ## sample size in simulation-based approach
      set.seed(seed) # set seed to make sure bootstrap replicate (draws) are correlated across functional domains
      for(i in 1:length(qn)){
        est_bs <- t(betaHat_boot[i,,])
        fit_fpca <- suppressWarnings( refund::fpca.face(est_bs, knots = nknots_fpca) ) # suppress sqrt(Eigen$values) NaNs
        ## extract estimated eigenfunctions/eigenvalues
        phi <- fit_fpca$efunctions
        lambda <- fit_fpca$evalues
        K <- length(fit_fpca$evalues)

        ## simulate random coefficients
        theta <- matrix(stats::rnorm(N*K), nrow=N, ncol=K) # generate independent standard normals
        if(K == 1){
          theta <- theta * sqrt(lambda) # scale to have appropriate variance
          X_new <- tcrossprod(theta, phi) # simulate new functions
        }else{
          theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
          X_new <- tcrossprod(theta, phi) # simulate new functions
        }
        x_sample <- X_new + t(fit_fpca$mu %o% rep(1,N)) # add back in the mean function
        Sigma_sd <- Rfast::colVars(x_sample, std = TRUE, na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
        x_mean <- colMeans(est_bs)
        un <- rep(NA, N)
        for(j in 1:N){
          un[j] <- max(abs((x_sample[j,] - x_mean)/Sigma_sd))
        }
        qn[i] <- stats::quantile(un, 0.95)
      }

      if(!silent)  message("Complete! \n -Use plot_fui() function to plot estimates \n -For more information, run the command:  ?plot_fui")

      return(list(betaHat = betaHat, betaHat.var = betaHat.var, qn = qn,
                  aic = AIC_mat, residuals = resids, bootstrap_samps = B,
                  argvals = argvals))

    }
  }else{
    # Only return betaHat when not calculating variance
    if(!silent)  message("Complete! \n -Use plot_fui() function to plot estimates \n -For more information, run the command:  ?plot_fui")

    return(list(betaHat = betaHat, argvals = argvals, aic = AIC_mat))
  }
}
