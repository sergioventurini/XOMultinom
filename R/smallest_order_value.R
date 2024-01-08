#' Survival function of the smallest order statistic for an equiprobable
#' multinomial distribution.
#'
#' This function calculates the survival function (i.e. the probability to be
#' greater than or equal to a given value) for the smallest order statistic
#' (i.e. the minimum) under an equiprobable multinomial distribution
#' assumption.
#'
#' @param t A length-one numeric vector indicating the value to compute the
#'   survival function for.
#' @param n A length-one integer vector indicating the number of independent
#'   balls.
#' @param m A length-one integer vector indicating the number of independent
#'   urns/cells.
#' @return A length-one numeric vector representing the probability of the
#'   smallest order statistic.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{highest_order_statistics}} for computing the
#'   CDF of the sum of the first \eqn{J} largest order statistics.
#' @seealso \code{\link{smallest_order_value}} for computing the CDF
#'   of the smallest order statistic.
#' @seealso \code{\link{range_probability}} for computing the CDF
#'   of the range.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' smallest_order_value(1, 10, 5) # P(N_(1) <= 1; n = 10, m = 5)
#'
#' @export
smallest_order_value <- function(t, n, m) {
  if (t > floor(n/m)) {
    P <- 0
    return(P)
  }
  if (t == 0) {
    P <- 1
    return(P)
  }
  P <- max_for_min(n, n, m, t)

  return(P)
}
