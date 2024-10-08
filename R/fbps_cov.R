# altered fbps() function from refund to increase GCV speed
# lambda1 == lambda2 for cov matrices b/c they're symmetric
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

  # return a smoothed matrix using fbps

  # data: a matrix
  # covariates: the list of data points for each dimension
  # knots: to specify either the number/numbers of  knots  or the
  # vector/vectors of knots for each dimension; defaults to 35
  # p: the degrees of B-splines
  # m: the order of difference penalty
  # lambda: the user-selected smoothing parameters
  # lscv: for leave-one-subject-out cross validation, the columns are
  # subjects
  # method: see "optim"
  # lower, upper, control: see "optim"

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
  #######################################################################
  Y = data

  ###################   precalculation for fbps smoothing  ##############

  List <- pspline_setting(x, xknots, p1, m1, periodicity[1])
  A1 <- List$A
  B1 <- List$B
  Bt1 <- Matrix(t(as.matrix(B1)))
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

  #################select optimal penalty ################################

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

  fbps_est =function(x) {

    lambda <- exp(x)
    ## two lambda's are the same
    if (length(lambda)==1)
    {
      lambda1 <- lambda
      lambda2 <- lambda
    }
    ## two lambda's are different
    if ( length(lambda) ==2) {
      lambda1 <- lambda[1]
      lambda2 <- lambda[2]
    }

    sigma2 <- 1/(1+lambda2*s2)
    sigma1 <- 1/(1+lambda1*s1)
    sigma2_sum <- sum(sigma2)
    sigma1_sum <- sum(sigma1)
    sigma <- Matrix::kronecker(sigma2,sigma1)
    sigma.2 <- Matrix::kronecker(sqrt(sigma2),sqrt(sigma1))

    Theta <- A01%*%diag(sigma1)%*%Ytilde
    Theta <- as.matrix(Theta%*%diag(sigma2)%*%t(A02))
    Yhat <- as.matrix(as.matrix(B1%*%Theta)%*%Bt2)

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

    if (search.grid ==T) {
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
