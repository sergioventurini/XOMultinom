#' PMF of the multinomial minimum at specified points
#'
#' Computes \eqn{P(\min(N_1, \ldots, N_m) = x)} at each element of \code{x}
#' for a multinomial random vector with \code{size} trials and cell
#' probabilities \code{prob}.  Returns a plain numeric vector, following the
#' same conventions as \code{\link[stats]{dbinom}} and
#' \code{\link[stats]{dnorm}}.
#'
#' @param x Integer vector of values at which to evaluate the PMF.
#' @param size Integer number of trials \eqn{n}.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalised to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param log.p Logical; if \code{TRUE}, log-probabilities are returned.
#'   Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation. Defaults to \code{TRUE}.
#'
#' @details
#' The function first checks whether \code{prob} corresponds to the
#' equiprobable case and then applies either the Bonetti et al.\ (2019)
#' algorithm or the Corrado (2011) algorithm accordingly.
#'
#' For the full distribution object (suitable for plotting, summaries, or
#' repeated evaluation), use \code{\link{minmultinomcdf}} directly.
#'
#' @return A numeric vector of the same length as \code{x}, containing
#'   \eqn{P(\min(N_1, \ldots, N_m) = x)} (or log-probabilities if
#'   \code{log.p = TRUE}).  Points outside the support \eqn{\{0, \ldots, n\}}
#'   return 0 (or \code{-Inf} on the log scale).
#'
#' @examples
#' m <- 4
#' n <- 60
#' probs <- rep(1 / m, m)
#'
#' # Evaluate at specific points -- plain numeric output, like dbinom()
#' dminmultinom(x = c(10, 12, 15), size = n, prob = probs)
#'
#' # Log scale
#' dminmultinom(x = c(10, 12, 15), size = n, prob = probs, log.p = TRUE)
#'
#' # For the full distribution object use minmultinomcdf():
#' Fmin <- minmultinomcdf(size = n, prob = probs)
#' plot(Fmin)
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
#' \code{\link{minmultinomcdf}} for the full distribution object,
#' \code{\link{pminmultinom}} for the CDF at specific points,
#' \code{\link{dmaxmultinom}} for the PMF of the maximum, and
#' \code{\link{drangemultinom}} for the PMF of the range.
#'
#' @export
dminmultinom <- function(x, size, prob, log.p = FALSE, verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob / s
  i0   <- prob == 0
  if (any(i0))
    prob <- prob[!i0]

  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- dminmultinom_bonetti(x, size, prob, FALSE, verbose)
  } else {
    res <- dminmultinom_corrado(x, size, prob, FALSE, verbose)
  }

  if (any(res < 0)) res[res < 0] <- 0
  if (any(res > 1)) res[res > 1] <- 1

  if (log.p) res <- log(res)

  res
}
