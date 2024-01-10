#' Utility function for computing recursively a sum.
#'
#' This is an auxiliary function to compute a sum recursively.
#'
#' @param t A length-one numeric vector indicating the value to compute the
#'   function for.
#' @param n A length-one integer vector indicating the number of independent
#'   balls.
#' @param m A length-one integer vector indicating the number of independent
#'   urns/cells.
#' @param J A length-one integer vector indicating the number of largest order
#'   statistics.
#' @param sum_depth A length-one integer vector indicating the summation depth.
#' @param cur_depth A length-one integer vector indicating the current
#'   summation depth.
#' @param rangeArg Integer vector indicating the range of the summation.
#' @return A length-one numeric vector representing the summation value.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{highest_order_statistics}} for computing the CDF
#'   of the sum for the first \eqn{J} largest order statistics for an
#'   equiprobable multinomial distribution.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' highest_order_statistics(6, 10, 5, 2)
#'
#' @export
recursive_sum <- function(t, n, m, J, sum_depth, cur_depth, rangeArg) {
  S <- 0
  cur_range <- numeric(2)

  if (cur_depth <= sum_depth) { # either increment summation depth or calculate the term
    if (cur_depth == 1) {
      cur_range[1] <- floor(t/J + 1)
      cur_range[2] <- t - sum_depth + 1
    } else {
      cur_range[1] <- floor((t - sum(rangeArg[1:(cur_depth - 1)]))/(J - cur_depth + 1) + 1)
      cur_range[2] <- min(c(rangeArg[cur_depth - 1], t - sum(rangeArg[1:(cur_depth - 1)])))
    }
    for (r in cur_range[1]:cur_range[2]) {
      rangeArg[cur_depth] <- r
      S <- S + recursive_sum(t, n, m, J, sum_depth, cur_depth + 1, rangeArg)
    }
  } else {
    prob_arg <- floor((t - sum(rangeArg))/(J - sum_depth))
    temp_p <- max_order_statistic(prob_arg, n - sum(rangeArg), m - sum_depth)
    common_term <- lgamma(n + 1) + lgamma(m + 1) - n*log(m)
    coef <- (n - sum(rangeArg))*log(m - sum_depth) - lgamma(m - sum_depth + 1) - lgamma(n - sum(rangeArg) + 1)
    for (k in 1:length(rangeArg)) {
      coef <- coef - lgamma(rangeArg[k] + 1)
    }
    equal_statistics <- unique(rangeArg)
    for (k in 1:length(equal_statistics)) {
      temp <- length(which(rangeArg == equal_statistics[k]))
      coef <- coef - lgamma(temp + 1)
    }
    S <- temp_p*exp(common_term + coef)
  }

  return(S)
}
