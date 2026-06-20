#' PMF of the sum of the \eqn{J} largest multinomial order statistics at
#' specified points
#'
#' Computes \eqn{P(S_J = x)}, where
#' \eqn{S_J = \sum_{j=1}^J N_{\langle j \rangle}}, at each element of \code{x}
#' for a multinomial random vector with \code{size} trials and equal cell
#' probabilities \code{prob}.  Returns a plain numeric vector, following the
#' same conventions as \code{\link[stats]{dbinom}} and
#' \code{\link[stats]{dnorm}}.
#'
#' @param x Integer vector of values at which to evaluate the PMF.
#' @param size Integer number of trials \eqn{n}.
#' @param prob Numeric vector of non-negative \emph{equal} cell probabilities
#'   (only the equiprobable case is implemented). Values are internally
#'   normalised to sum to 1.
#' @param J Integer; number of largest order statistics to sum. Defaults to
#'   \code{2}.
#' @param log.p Logical; if \code{TRUE}, log-probabilities are returned.
#'   Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation. Defaults to \code{TRUE}.
#'
#' @details
#' Only the equiprobable case (\code{prob} proportional to a constant vector)
#' is currently supported.
#'
#' For the full distribution object (suitable for plotting, summaries, or
#' repeated evaluation), use \code{\link{Jlargemultinomcdf}} directly.
#'
#' @return A numeric vector of the same length as \code{x}, containing
#'   \eqn{P(S_J = x)} (or log-probabilities if \code{log.p = TRUE}).  Points
#'   outside the support \eqn{\{0, \ldots, n\}} return 0 (or \code{-Inf} on
#'   the log scale).
#'
#' @examples
#' m <- 4
#' n <- 60
#' probs <- rep(1 / m, m)
#' J <- 3
#'
#' # Evaluate at specific points -- plain numeric output, like dbinom()
#' dJlargemultinom(x = c(30, 35, 40), size = n, prob = probs, J = J)
#'
#' # Log scale
#' dJlargemultinom(x = c(30, 35, 40), size = n, prob = probs, J = J,
#'                 log.p = TRUE)
#'
#' # For the full distribution object use Jlargemultinomcdf():
#' FJ <- Jlargemultinomcdf(size = n, prob = probs, J = J)
#' plot(FJ)
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' @seealso
#' \code{\link{Jlargemultinomcdf}} for the full distribution object,
#' \code{\link{pJlargemultinom}} for the CDF at specific points,
#' \code{\link{dmaxmultinom}} for the PMF of the maximum, and
#' \code{\link{dminmultinom}} for the PMF of the minimum.
#'
#' @export
dJlargemultinom <- function(x, size, prob, J = 2, log.p = FALSE,
                            verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob / s
  i0   <- prob == 0
  if (any(i0))
    prob <- prob[!i0]

  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- dJlargemultinom_bonetti(x, size, prob, J, FALSE, verbose)
  } else {
    stop("not available for the non-equiprobable case.")
  }

  if (any(res < 0)) res[res < 0] <- 0
  if (any(res > 1)) res[res > 1] <- 1

  if (log.p) res <- log(res)

  res
}
