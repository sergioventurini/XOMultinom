#' CDF of the multinomial maximum at specified points
#'
#' Computes the cumulative distribution function of the maximum cell count of a
#' multinomial random vector with arbitrary cell probabilities.
#'
#' @param x Numeric vector of values at which to evaluate the CDF.
#' @param size Integer number of trials \eqn{n}.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalised to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param lower.tail Logical; if \code{TRUE} (default),
#'   \eqn{P(\max(N_1, \ldots, N_m) \le x)} is returned; otherwise
#'   \eqn{P(\max(N_1, \ldots, N_m) > x)}.
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
#'   \eqn{P(\max(N_1, \ldots, N_m) \le x)} (or its complement or log, according
#'   to \code{lower.tail} and \code{log.p}). Values outside the support are
#'   handled consistently with base R: \code{x < 0} gives 0 and \code{x > n}
#'   gives 1 (before \code{lower.tail}/\code{log.p} transformations).
#'
#' @examples
#' m <- 4
#' n <- 60
#' probs <- rep(1 / m, m)
#' xseq <- 0:n
#'
#' cdfmax <- pmaxmultinom(x = xseq, size = n, prob = probs)
#' cdfmax
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
#' \code{\link{maxmultinomcdf}} for the full distribution object,
#' \code{\link{pminmultinom}} for the CDF of the minimum,
#' \code{\link{dmaxmultinom}} for the PMF of the maximum, and
#' \code{\link{dminmultinom}} for the PMF of the minimum.
#'
#' @export
pmaxmultinom <- function(x, size, prob, lower.tail = TRUE, log.p = FALSE,
                         verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob/s
  i0 <- prob == 0
  if (any(i0))
    prob <- prob[!i0]
  
  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- pmaxmultinom_bonetti(x, size, prob, FALSE, verbose)
  } else {
    res <- pmaxmultinom_corrado(x, size, prob, FALSE, verbose)
  }

  if (any(res < 0)) res[res < 0] <- 0
  if (any(res > 1)) res[res > 1] <- 1

  if (!lower.tail) res <- 1 - res
  if (log.p)       res <- log(res)

  res
}

#' Distribution object for the multinomial maximum count
#'
#' Constructs an \code{xomultinom_dist} object containing the exact CDF
#' of the maximum cell count \eqn{\max(N_1, \ldots, N_m)} of a multinomial
#' random vector, evaluated over its full support \eqn{\{0, 1, \ldots, n\}}.
#' The returned object can be passed to \code{plot()}, \code{autoplot()},
#' \code{summary()}, and \code{as.data.frame()}, and its CDF and PMF values can
#' be extracted with \code{pmaxmultinom()} and \code{dmaxmultinom()}.
#'
#' @param size Integer number of trials \eqn{n}.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalised to sum to 1. Categories with zero probability are
#'   removed before computation.
#' @param verbose Logical; if \code{TRUE}, displays progress information during
#'   the computation. Defaults to \code{TRUE}.
#'
#' @details
#' \code{maxmultinomcdf()} is the \emph{distribution constructor}: it fixes
#' \code{size} and \code{prob}, performs the (potentially expensive) exact
#' computation once over the full support, and returns a self-contained
#' \code{xomultinom_dist} object.  The companion functions
#' \code{\link{pmaxmultinom}} and \code{\link{dmaxmultinom}} provide the CDF or
#' PMF values at the requested points \code{x}, returning a plain numeric
#' vector in the same style as \code{\link[stats]{pnorm}} and
#' \code{\link[stats]{dnorm}}.
#'
#' Use \code{maxmultinomcdf()} when you need the full distribution object (e.g.,
#' for plotting or for evaluating the CDF at many points without repeating the
#' underlying computation).  Use \code{\link{pmaxmultinom}} or
#' \code{\link{dmaxmultinom}} when you need a numeric vector at specific
#' quantiles, in the same way you would use \code{pnorm()} or \code{dnorm()}.
#'
#' The function dispatches automatically to the Bonetti et al. (2019) recursive
#' algorithm (equiprobable case) or the Corrado (2011) matrix algorithm
#' (general case).
#'
#' @return An object of class \code{xomultinom_dist} with components \code{x}
#'   (full integer support \eqn{0, \ldots, n}), \code{values} (CDF values),
#'   \code{stat = "max"}, \code{type = "cdf"}, \code{size}, \code{prob}, and
#'   \code{log = FALSE}.
#'
#' @examples
#' m <- 4; n <- 60
#' probs <- rep(1 / m, m)
#'
#' # Distribution constructor: compute once, reuse freely
#' Fmax <- maxmultinomcdf(size = n, prob = probs)
#' plot(Fmax)
#' summary(Fmax)
#'
#' # Standard p*/d* interface: plain numeric output
#' pmaxmultinom(x = c(18, 20, 22), size = n, prob = probs)
#' dmaxmultinom(x = c(18, 20, 22), size = n, prob = probs)
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
#' \code{\link{pmaxmultinom}} for the CDF at specific points (numeric output),
#' \code{\link{dmaxmultinom}} for the PMF at specific points (numeric output),
#' \code{\link{minmultinomcdf}} and \code{\link{rangemultinomcdf}} for the analogous
#' constructors for the minimum and the range.
#'
#' @export
maxmultinomcdf <- function(size, prob, verbose = TRUE) {
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob / s
  i0   <- prob == 0
  if (any(i0))
    prob <- prob[!i0]

  xseq <- 0:size

  if (length(prob) > 0 && length(unique(prob)) == 1) {
    res <- pmaxmultinom_bonetti(xseq, size, prob, FALSE, verbose)
  } else {
    res <- pmaxmultinom_corrado(xseq, size, prob, FALSE, verbose)
  }

  if (any(res < 0)) res[res < 0] <- 0
  if (any(res > 1)) res[res > 1] <- 1

  new_xomultinom_dist(x      = xseq,
                      values = res,
                      stat   = "max",
                      type   = "cdf",
                      size   = size,
                      prob   = prob,
                      log    = FALSE)
}
