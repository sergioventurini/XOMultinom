#' CDF of the multinomial range at specified points
#'
#' Computes the cumulative distribution function of the range
#' \eqn{R = \max(N_1, \ldots, N_m) - \min(N_1, \ldots, N_m)}
#' for a multinomial random vector with arbitrary cell probabilities.
#'
#' @param x Numeric vector of values at which to evaluate the CDF.
#' @param size Integer number of trials \eqn{n}.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalised to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param lower.tail Logical; if \code{TRUE} (default),
#'   \eqn{P(R \le x)} is returned; otherwise \eqn{P(R > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned on the log
#'   scale. Defaults to \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation. Defaults to \code{TRUE}.
#'
#' @details The function first checks whether `prob` corresponds to the
#'   equiprobable case and then applies either the Bonetti et al. (2019)
#'   algorithm or the Corrado (2011) algorithm accordingly.
#'
#' @return A numeric vector of the same length as \code{x}, containing
#'   \eqn{P(R \le x)} (or its complement or log, according to
#'   \code{lower.tail} and \code{log.p}). Values outside the support are
#'   handled consistently with base R: \code{x < 0} gives 0 and \code{x > n}
#'   gives 1 (before \code{lower.tail}/\code{log.p} transformations).
#'
#' @examples
#' m <- 4
#' n <- 60
#' probs <- rep(1 / m, m)
#' xseq <- 0:n
#'
#' cdfrange <- prangemultinom(x = xseq, size = n, prob = probs)
#' cdfrange
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
#' \code{\link{prangemultinom}} for the CDF at specific points (numeric output),
#' \code{\link{drangemultinom}} for the PMF at specific points (numeric output),
#' \code{\link{maxmultinomcdf}} and \code{\link{minmultinomcdf}} for the analogous
#' constructors for the maximum and the minimum.
#'
#' @export
prangemultinom <- function(x, size, prob, lower.tail = TRUE, log.p = FALSE,
                           verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob/s
  i0 <- prob == 0
  if (any(i0))
    prob <- prob[!i0]

  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- prangemultinom_bonetti(x, size, prob, FALSE, verbose)
  } else {
    res <- prangemultinom_corrado(x, size, prob, FALSE, verbose)
  }

  if (any(res < 0)) res[res < 0] <- 0
  if (any(res > 1)) res[res > 1] <- 1

  if (!lower.tail) res <- 1 - res
  if (log.p)       res <- log(res)

  res
}

#' Distribution object for the multinomial range
#'
#' Constructs an \code{xomultinom_dist} object containing the exact PMF and CDF
#' of the range \eqn{R = \max(N_1, \ldots, N_m) - \min(N_1, \ldots, N_m)} of a
#' multinomial random vector, evaluated over its full support
#' \eqn{\{0, 1, \ldots, n\}}.  The returned object can be passed to
#' \code{plot()}, \code{autoplot()}, \code{summary()}, and
#' \code{as.data.frame()}, and its CDF and PMF values can be extracted with
#' \code{prangemultinom()} and \code{drangemultinom()}.
#'
#' @param size Integer number of trials \eqn{n}.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalised to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation. Defaults to \code{TRUE}.
#'
#' @details
#' \code{rangemultinomcdf()} is the \emph{distribution constructor}: it fixes
#' \code{size} and \code{prob}, performs the exact computation once over the
#' full support, and returns a self-contained \code{xomultinom_dist} object.
#' The companion functions \code{\link{prangemultinom}} and
#' \code{\link{drangemultinom}} are lightweight wrappers that call
#' \code{rangemultinomcdf()} internally and extract the CDF or PMF values at the
#' requested points \code{x}, returning a plain numeric vector in the same
#' style as \code{\link[stats]{pnorm}} and \code{\link[stats]{dnorm}}.
#'
#' Use \code{rangemultinomcdf()} when you need the full distribution object (e.g.,
#' for plotting or for evaluating the CDF at many points without repeating the
#' underlying computation).  Use \code{\link{prangemultinom}} or
#' \code{\link{drangemultinom}} when you need a numeric vector at specific
#' quantiles.
#'
#' The function dispatches automatically to the Bonetti et al. (2019) recursive
#' algorithm (equiprobable case) or the Corrado (2011) matrix algorithm
#' (general case).
#'
#' @return An object of class \code{xomultinom_dist} with components \code{x}
#'   (full integer support \eqn{0, \ldots, n}), \code{values} (CDF values),
#'   \code{stat = "range"}, \code{type = "cdf"}, \code{size}, \code{prob}, and
#'   \code{log = FALSE}.
#'
#' @examples
#' m <- 4; n <- 60
#' probs <- rep(1 / m, m)
#'
#' # Distribution constructor: compute once, reuse freely
#' Frange <- rangemultinomcdf(size = n, prob = probs)
#' plot(Frange)
#' summary(Frange)
#'
#' # Standard p*/d* interface: plain numeric output
#' prangemultinom(x = c(5, 10, 15), size = n, prob = probs)
#' drangemultinom(x = c(5, 10, 15), size = n, prob = probs)
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
#' \code{\link{prangemultinom}} for the CDF at specific points (numeric output),
#' \code{\link{drangemultinom}} for the PMF at specific points (numeric output),
#' \code{\link{maxmultinomcdf}} and \code{\link{minmultinomcdf}} for the analogous
#' constructors for the maximum and the minimum.
#'
#' @export
rangemultinomcdf <- function(size, prob, verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob / s
  i0   <- prob == 0
  if (any(i0))
    prob <- prob[!i0]
 
  xseq <- 0:size
 
  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- prangemultinom_bonetti(xseq, size, prob, FALSE, verbose)
  } else {
    res <- prangemultinom_corrado(xseq, size, prob, FALSE, verbose)
  }
 
  if (any(res < 0)) res[res < 0] <- 0
  if (any(res > 1)) res[res > 1] <- 1
 
  new_xomultinom_dist(x      = xseq,
                      values = res,
                      stat   = "range",
                      type   = "cdf",
                      size   = size,
                      prob   = prob,
                      log    = FALSE)
}
