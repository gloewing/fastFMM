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
#'
#' @export

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
