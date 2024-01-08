#' Survival function of smallest order statistic for a multinomial
#' distribution.
#'
#' This function calculates the survival function (i.e. the probability to be
#' greater than or equal to a given value) for the smallest order statistic
#' (i.e. the minimum) under an equiprobable multinomial distribution
#' assumption.
#'
#' @param t_max A length-one numeric vector indicating the value to compute the
#'   survival function for.
#' @param n A length-one integer vector indicating the number of independent
#'   balls.
#' @param m A length-one integer vector indicating the number of independent
#'   urns/cells.
#' @param t A length-one numeric vector indicating the value to compute the
#'   survival function for.
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
#' smallest_order_value(3, 10, 5) # P(N_(1) <= 3; n = 10, m = 5)
#'
#' @export
max_for_range <- function(t_max, n, m, prev, t) {
  aux <- 0
  if (identical(c(t_max, n, m), prev)) {
    aux <- 0
    return(aux)
  }
  if (prev[1] + 1 - t_max > t) {
    if (n == 0 & m == 0) {
      aux <- 1
      return(aux)
    } else {
      aux <- 0
      return(aux)
    }
  } else {
    if (n == 0 & m == 0) {
      aux <- 1
      return(aux)
    }
    if (t_max == 1) {
      if (m == n) {
        aux <- exp(lgamma(m + 1) - lgamma(m - n + 1) - n*log(m)) # explicit calculation of P(n<1> <= 1)
        return(aux)
      } else {
        aux <- 0
        return(aux)
      }
    }
    if (n == 0 & m != 0) {
      aux <- 0
      return(aux)
    }
    common_term <- lgamma(n + 1) + lgamma(m + 1) - n*log(m)
    LowSum <- max(0, n - t_max*m + m)
    UpSum <- floor(n/t_max)
    if (LowSum <= UpSum) {
      for (q in LowSum:UpSum) {
        summ_term <- -q*lgamma(t_max + 1) - lgamma(q + 1) - lgamma(m - q + 1) - lgamma(n - t_max*q + 1)
        if (m == q) {
          summ_term_nominator <- 0
        } else {
          summ_term_nominator <- (n - t_max*q)*log(m - q)
        }
        coef <- exp(common_term + summ_term + summ_term_nominator)
        temp <- max_for_range(t_max - 1, m - t_max*q, m - q, prev, t)
        aux <- aux + coef*temp
      }
    }
  }

  return(aux)
}
