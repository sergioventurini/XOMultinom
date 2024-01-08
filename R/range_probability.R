#' CDF of the range for an equiprobable multinomial distribution.
#'
#' This function calculates the cumulative distribution function (CDF)
#' of the range under the assumption of an equiprobable multinomial
#' distribution.
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
#' @seealso \code{\link{highest_order_statistics}} for computing the sum of
#'   the first \eqn{J} largest order statistics.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' range_probability(1, 7, 10)
#'
#' @export
range_probability <- function(t, n, m) {
  P <- 0
  t <- floor(t)
  if (t > n) {
    P <- 1
    return(P)
  }
  P <- max_order_statistic(t, n, m)
  prev <- c(t, n, m)
  for (t_max in (t + 1):n) {
    aux <- max_for_range(t_max, n, m, prev, t)
    P <- P + aux
    prev <- c(t_max, n, m)
  }

  return(P)
}
