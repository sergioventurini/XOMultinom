#' CDF of the maximum for a multinomial distribution
#'
#' Computes the cumulative distribution function of the maximum cell count of a
#' multinomial random vector with arbitrary cell probabilities.
#'
#' @param x Numeric vector of values at which to evaluate the CDF.
#' @param size Integer number of trials.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalized to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param log Logical; if \code{TRUE}, returns log-probabilities.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation.
#'
#' @return Numeric vector of the same length as \code{x} containing
#'   \eqn{P(\max(X_1, \ldots, X_K) \le x)} values, or log-probabilities if
#'   \code{log = TRUE}.
#'
#' @examples
#' k <- 4
#' n <- 60
#' probs <- rep(1 / k, k)
#' xseq <- 0:n
#'
#' cdfmax <- pmaxmultinom(x = xseq, size = n, prob = probs)
#' plot(cdfmax)
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' @seealso
#' \code{\link{pminmultinom}} for the CDF of the minimum,
#' \code{\link{dmaxmultinom}} for the PMF of the maximum, and
#' \code{\link{dminmultinom}} for the PMF of the minimum.
#'
#' @export
pmaxmultinom <- function(x, size, prob, log = FALSE, verbose = FALSE) {
  k <- length(prob)
  xlen <- length(x)
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob/s
  # x <- as.integer(x + 0.5)
  i0 <- prob == 0
  if (any(i0)) {
    prob <- prob[!i0]
    k <- length(prob)
  }
  
  res <- pmaxmultinom_corrado(x, size, prob, log, verbose)

  if (any(res < 0))
    res[res < 0] <- 0
  if (any(res > 1))
    res[res > 1] <- 1

  return(new_xomultinom_dist(x      = x,
                             values = res,
                             stat   = "max",
                             type   = "cdf",
                             size   = size,
                             prob   = prob,
                             log    = log))
}
