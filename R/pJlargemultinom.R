#' CDF of the sum of \eqn{J} largest order statistics for a multinomial
#' distribution evaluated at specified points
#'
#' Computes the cumulative distribution function of the sum of \eqn{J}
#' largest order statistics, \eqn{S_J = \sum_{j=1}^J N_{\langle j\rangle}},
#' for a multinomial random vector with equal cell probabilities.
#'
#' @param x Numeric vector of values at which to evaluate the CDF.
#' @param size Integer number of trials \eqn{n}.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalised to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param J Integer; number of largest order statistics to sum. Defaults to
#'   \code{2}.
#' @param lower.tail Logical; if \code{TRUE} (default),
#'   \eqn{P(S_J \le x)} is returned; otherwise \eqn{P(S_J > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned on the log
#'   scale. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation. Defaults to \code{TRUE}.
#'
#' @details The function only implements the equiprobable case.
#'
#' @return A numeric vector of the same length as \code{x}, containing
#'   \eqn{P(S_J \le x)} (or its complement or log, according to
#'   \code{lower.tail} and \code{log.p}). Values outside the support are
#'   handled consistently with base R: \code{x < 0} gives 0 and \code{x > n}
#'   gives 1 (before \code{lower.tail}/\code{log.p} transformations).
#'
#' @examples
#' m <- 4
#' n <- 60
#' probs <- rep(1 / m, m)
#' J <- 3
#' xseq <- 0:n
#'
#' cdflarge <- pJlargemultinom(x = xseq, size = n, prob = probs, J = J)
#' cdflarge
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' @seealso
#' \code{\link{Jlargemultinomcdf}} for the full distribution object,
#' \code{\link{dJlargemultinom}} for the PMF,
#' \code{\link{qJlargemultinom}} for quantiles, and
#' \code{\link{rJlargemultinom}} for random generation.
#'
#' @export
pJlargemultinom <- function(x, size, prob, J = 2, lower.tail = TRUE, log.p = FALSE,
                            verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob/s
  i0 <- prob == 0
  if (any(i0))
    prob <- prob[!i0]
  
  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- pJlargemultinom_bonetti(x, size, prob, J, FALSE, verbose)
  } else {
    stop("not available for the non-equiprobable case.")
  }

  if (any(res < 0)) res[res < 0] <- 0
  if (any(res > 1)) res[res > 1] <- 1

  if (!lower.tail) res <- 1 - res
  if (log.p)       res <- log(res)

  res
}

#' Distribution object for the sum of the \eqn{J} largest multinomial order
#' statistics
#'
#' Constructs an \code{xomultinom_dist} object containing the exact PMF and CDF
#' of \eqn{S_J = \sum_{j=1}^J N_{\langle j \rangle}}, the sum of the \eqn{J}
#' largest order statistics of a multinomial random vector, evaluated over its
#' full support \eqn{\{0, 1, \ldots, n\}}.  The returned object can be passed
#' to \code{plot()}, \code{autoplot()}, \code{summary()}, and
#' \code{as.data.frame()}, and its CDF and PMF values can be extracted with
#' \code{pJlargemultinom()} and \code{dJlargemultinom()}.
#'
#' @param size Integer number of trials \eqn{n}.
#' @param prob Numeric vector of non-negative \emph{equal} cell probabilities
#'   (only the equiprobable case is implemented). Values are internally
#'   normalised to sum to 1.
#' @param J Integer; number of largest order statistics to sum. Defaults to
#'   \code{2}.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation. Defaults to \code{TRUE}.
#'
#' @details
#' \code{Jlargemultinomcdf()} is the \emph{distribution constructor}: it fixes
#' \code{size}, \code{prob}, and \code{J}, performs the exact computation once
#' over the full support, and returns a self-contained \code{xomultinom_dist}
#' object.  The companion functions \code{\link{pJlargemultinom}} and
#' \code{\link{dJlargemultinom}} are lightweight wrappers that call
#' \code{Jlargemultinomcdf()} internally and extract the CDF or PMF values at the
#' requested points \code{x}, returning a plain numeric vector in the same
#' style as \code{\link[stats]{pnorm}} and \code{\link[stats]{dnorm}}.
#'
#' Only the equiprobable case (\code{prob} proportional to a constant vector)
#' is currently supported.
#'
#' @return An object of class \code{xomultinom_dist} with components \code{x}
#'   (full integer support \eqn{0, \ldots, n}), \code{values} (CDF values),
#'   \code{stat = "J_largest"}, \code{type = "cdf"}, \code{size}, \code{prob},
#'   and \code{log = FALSE}.
#'
#' @examples
#' m <- 4; n <- 60; J <- 3
#' probs <- rep(1 / m, m)
#'
#' # Distribution constructor: compute once, reuse freely
#' FJ <- Jlargemultinomcdf(size = n, prob = probs, J = J)
#' plot(FJ)
#' summary(FJ)
#'
#' # Standard p*/d* interface: plain numeric output
#' pJlargemultinom(x = c(30, 35, 40), size = n, prob = probs, J = J)
#' dJlargemultinom(x = c(30, 35, 40), size = n, prob = probs, J = J)
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' @seealso
#' \code{\link{pJlargemultinom}} for the CDF at specific points (numeric
#' output), \code{\link{dJlargemultinom}} for the PMF at specific points,
#' \code{\link{maxmultinomcdf}}, \code{\link{minmultinomcdf}}, and
#' \code{\link{rangemultinomcdf}} for the analogous constructors.
#'
#' @export
Jlargemultinomcdf <- function(size, prob, J = 2, verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob / s
  i0   <- prob == 0
  if (any(i0))
    prob <- prob[!i0]
 
  if (!(length(prob) > 0 && length(unique(prob)) == 1))
    stop("not available for the non-equiprobable case.")
 
  xseq <- 0:size
  res  <- pJlargemultinom_bonetti(xseq, size, prob, J, FALSE, verbose)
 
  if (any(res < 0)) res[res < 0] <- 0
  if (any(res > 1)) res[res > 1] <- 1
 
  new_xomultinom_dist(x      = xseq,
                      values = res,
                      stat   = "J_largest",
                      type   = "cdf",
                      size   = size,
                      prob   = prob,
                      log    = FALSE)
}
