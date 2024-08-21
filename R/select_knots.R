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
#' @export

# TODO: figure out what `t` is and whether p is num or int

select_knots <- function(
  t,
  knots = 10,
  p = 3,
  option = "equally-spaced"){


  qs <- seq(0, 1, length = knots+1)

  if(option == "equally-spaced"){
    knots <- (max(t)-min(t))*qs + min(t)
  }
  if(option == "quantile") {
    knots <- as.vector(stats::quantile(t,qs))
  }

  K <- length(knots)
  knots_left <- 2*knots[1] - knots[p:1+1] # AX: Why isn't this p:2
  knots_right <- 2*knots[K] - knots[K - (1:p)]

  return(c(knots_left, knots, knots_right))
}
