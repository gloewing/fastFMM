      G_estimate <- function(data, out_index, data_cov, ztlist, designmat,
                             betaHat, HHat, RE_table, non_neg = 1, MoM = 2,
                             silent = TRUE){

        if(silent == FALSE) print("Step 3.1.1: Method of Moments Covariance Estimator")

        GTilde <- array(NA, dim = c(nrow(HHat), L, L))
        idx_lst <- data_cov$idx_lst # indices of beta to average over
        z_names <- names(ztlist)
        Z <- as.matrix(data_cov$Z_orig) # this has concatenated design matrix and sums over columns for ID variables #do.call(cbind, ztlist)

        # first part of OLS
        if(MoM == 2){
          XTXX <- as.matrix( tcrossprod( MASS::ginv( as.matrix(crossprod(data_cov$Z)) ), data_cov$Z) ) # first part of OLS expression
        }else if(MoM == 1){
          # function to join matrices
          mat_concat <- function(yy, xx){
            if(length(xx) > 1){
              return( rowSums(yy[,xx]) )
            }else{
              return( yy[,xx] )
            }
          }
          ZZ <- do.call( cbind, lapply( idx_lst,  function(xx) mat_concat(yy=data_cov$Z, xx=xx) ) ) # sum across columns of Z associated with same random effect
          XTXX <- tcrossprod( MASS::ginv( crossprod(ZZ) ), ZZ) # first part of OLS expression
          rm(ZZ)
        }

        # save design matrix for random effects
        idx_vec_zlst <- sapply(ztlist, ncol) # number of cols in each matrix
        idx_vec_zlst <- c(0, cumsum(idx_vec_zlst)) # vector of indices

        d_temp <- data[,out_index]
        if(MoM == 2){
          for(i in 1:L){
            YYi <- (d_temp[,i] - designmat %*% betaHat[,i])
            for(j in i:L){
              YYj <- YYi * (d_temp[,j] - designmat %*% betaHat[,j]) # outcome of product of residuals
              bHat <- XTXX %*% YYj  # coefficients from OLS with pseudo-inverse
              GTilde[,i,j] <- GTilde[,j,i] <- sapply( idx_lst,  function(x) mean(bHat[x]) )
            }
          }
          rm(bHat)
        }else if(MoM == 1){
          for(i in 1:L){
            YYi <- (d_temp[,i] - designmat %*% betaHat[,i])
            for(j in i:L){
              YYj <- YYi * (d_temp[,j] - designmat %*% betaHat[,j]) # outcome of product of residuals
              GTilde[,i,j] <- GTilde[,j,i] <- XTXX %*% YYj # coefficients from OLS with pseudo-inverse
            }
          }
        }
        rm(YYi, YYj, XTXX)

        # non-negative least squares for estimation of non-diagonal terms
        if(non_neg != 0){

          if(MoM == 1 & non_neg == 2){
            message("Method of Moments approach 1 estimator can only use NNLS estimation scheme 1. Proceeding with NNLS-1")
            non_neg <- 1
          }

          GTilde <- cov.nnls(data_cov = data_cov,
                             RE_table = RE_table,
                             idx_lst = idx_lst,
                             d_temp = d_temp,
                             designmat = designmat,
                             betaHat = betaHat,
                             GTilde = GTilde,
                             non_neg = non_neg,
                             silent = silent)
        }

        dimnames(GTilde)[[1]] <- rownames(HHat) # use names so lfosr_cov_organize() function below knows how to organize sub matrices based on names

        return(GTilde)
      }
