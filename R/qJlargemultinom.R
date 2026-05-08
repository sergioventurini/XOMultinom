#' Quantile function of the sum of \eqn{J} largest order statistics for a
#' multinomial distribution
#'
#' Computes exact quantiles of the distribution of the sum of the \eqn{J}
#' largest order statistics \eqn{S_J = \sum_{j=1}^{J} N_{\langle j \rangle}}
#' of a multinomial random vector with equal cell probabilities, by inverting
#' the exact CDF obtained from \code{\link{pJlargemultinom}}.
#'
#' @param p Numeric vector of probabilities (or log-probabilities if
#'   \code{log.p = TRUE}) at which to evaluate the quantile function.
#' @param size Integer number of trials.
#' @param prob Numeric vector of non-negative, equal cell probabilities.
#'   Only the equiprobable case is supported; a non-equiprobable \code{prob}
#'   will raise an error (propagated from \code{\link{pJlargemultinom}}).
#' @param J Integer number of largest order statistics to sum. Defaults to
#'   \code{2}.
#' @param lower.tail Logical; if \code{TRUE} (default),
#'   \eqn{Q(p) = \min\{x : F(x) \ge p\}}; if \code{FALSE},
#'   \eqn{Q(p) = \min\{x : F(x) \ge 1 - p\}}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is taken to be on the log
#'   scale. Defaults to \code{FALSE}.
#'
#' @return Integer vector of the same length as \code{p} containing the
#'   corresponding exact quantiles of \eqn{S_J}.
#'
#' @details
#'   The function obtains the exact CDF over the full support
#'   \eqn{\{0, 1, \ldots, n\}} via a single vectorised call to
#'   \code{\link{pJlargemultinom}}.  The quantile is then located as the
#'   smallest support point whose CDF value meets or exceeds \code{p}.
#'   Only the equiprobable case is supported, consistent with
#'   \code{\link{pJlargemultinom}}.
#'
#' @examples
#' m <- 4
#' n <- 60
#' probs <- rep(1 / m, m)
#'
#' # Median and 95th percentile of S_3
#' qJlargemultinom(c(0.5, 0.95), size = n, prob = probs, J = 3)
#'
#' # Upper tail
#' qJlargemultinom(0.05, size = n, prob = probs, J = 3, lower.tail = FALSE)
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' @seealso
#' \code{\link{pJlargemultinom}} for the CDF,
#' \code{\link{dJlargemultinom}} for the PMF,
#' \code{\link{rJlargemultinom}} for random generation.
#'
#' @export
qJlargemultinom <- function(p, size, prob, J = 2,
                            lower.tail = TRUE, log.p = FALSE) {
  if (log.p)
    p <- exp(p)
  if (any(p < 0 | p > 1, na.rm = TRUE))
    stop("'p' must contain values in [0, 1].")
  if (!lower.tail)
    p <- 1 - p
 
  supp <- 0L:size
  cdf  <- pJlargemultinom(x = supp, size = size, prob = prob, J = J,
                          log = FALSE, verbose = FALSE)$values
 
  discrete_quantile(p, cdf, supp)
}
