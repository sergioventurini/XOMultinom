#' PMF of the sum of \eqn{J} largest order statistics for a multinomial
#' distribution
#'
#' Computes the probability mass function of the sum of \eqn{J}
#' largest order statistics for a multinomial random vector with
#' equal cell probabilities.
#'
#' @param x Numeric vector of values at which to evaluate the PMF.
#' @param size Integer number of trials.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalized to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param J Integer number of largest order statistics to consider. Defaults to
#'   \code{2}.
#' @param log Logical; if \code{TRUE}, returns log-probabilities.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation.
#'
#' @details The function only implements the equiprobable case.
#'
#' @return An object of class \code{xomultinom_dist} with fields \code{x},
#'   \code{values} (containing \eqn{P(S_J = x)},
#'   \eqn{S_J = \sum_{j=1}^J N_{\langle j\rangle}}, or log-probabilities if
#'   \code{log = TRUE}), \code{stat = "J_largest"}, \code{type = "pmf"},
#'   \code{size}, \code{prob}, and \code{log}.
#'
#' @examples
#' m <- 4
#' n <- 60
#' probs <- rep(1 / m, m)
#' J <- 3
#' xseq <- 0:n
#'
#' pmflarge <- dJlargemultinom(x = xseq, size = n, prob = probs, J = J)
#' plot(pmflarge)
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' @seealso
#' \code{\link{dmaxmultinom}} for the PMF of the maximum,
#' \code{\link{dminmultinom}} for the PMF of the minimum, and
#' \code{\link{drangemultinom}} for the PMF of the range.
#'
#' @export
dJlargemultinom <- function(x, size, prob, J = 2, log = FALSE, verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob/s
  i0 <- prob == 0
  if (any(i0))
    prob <- prob[!i0]
  
  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- dJlargemultinom_bonetti(x, size, prob, J, log, verbose)
  } else {
    stop("not available for the non-equiprobable case.")
  }

  if (!log) {
    if (any(res < 0)) res[res < 0] <- 0
    if (any(res > 1)) res[res > 1] <- 1
  } else {
    if (any(res > 0)) res[res > 0] <- 0
  }

  return(new_xomultinom_dist(x      = x,
                             values = res,
                             stat   = "J_largest",
                             type   = "pmf",
                             size   = size,
                             prob   = prob,
                             log    = log))
}
