#' PMF of the maximum for a multinomial distribution
#'
#' Computes the probability mass function of the maximum cell count of a
#' multinomial random vector with arbitrary cell probabilities.
#'
#' @param x Numeric vector of values at which to evaluate the PMF.
#' @param size Integer number of trials.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalized to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param log Logical; if \code{TRUE}, returns log-probabilities.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation.
#'
#' @details The function first checks whether `prob` corresponds to the
#'   equiprobable case and then applies either the Bonetti et al. (2019)
#'   algorithm or the Corrado (2011) algorithm accordingly.
#'
#' @return An object of class \code{xomultinom_dist} with fields \code{x},
#'   \code{values} (containing \eqn{P(\max(N_1, \ldots, N_m) = x)}, or
#'   log-probabilities if \code{log = TRUE}), \code{stat = "max"},
#'   \code{type = "pmf"}, \code{size}, \code{prob}, and \code{log}.
#'
#' @examples
#' m <- 4
#' n <- 60
#' probs <- rep(1 / m, m)
#' xseq <- 0:n
#'
#' pmfmax <- dmaxmultinom(x = xseq, size = n, prob = probs)
#' plot(pmfmax)
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' Corrado, C.J. (2011). The exact distribution of the maximum,
#' minimum and the range of Multinomial/Dirichlet and Multivariate
#' Hypergeometric frequencies. Statistical Computing, 21, 349--359.
#' \doi{10.1007/s11222-010-9174-3}
#'
#' @seealso
#' \code{\link{pmaxmultinom}} for the CDF of the maximum,
#' \code{\link{dminmultinom}} for the PMF of the minimum, and
#' \code{\link{drangemultinom}} for the PMF of the range.
#'
#' @export
dmaxmultinom <- function(x, size, prob, log = FALSE, verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob/s
  i0 <- prob == 0
  if (any(i0))
    prob <- prob[!i0]
  
  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- dmaxmultinom_bonetti(x, size, prob, log, verbose)
  } else {
    res <- dmaxmultinom_corrado(x, size, prob, log, verbose)
  }

  if (!log) {
    if (any(res < 0)) res[res < 0] <- 0
    if (any(res > 1)) res[res > 1] <- 1
  } else {
    if (any(res > 0)) res[res > 0] <- 0
  }

  return(new_xomultinom_dist(x      = x,
                             values = res,
                             stat   = "max",
                             type   = "pmf",
                             size   = size,
                             prob   = prob,
                             log    = log))
}
