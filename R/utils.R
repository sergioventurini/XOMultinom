#' Stable rounding function
#'
#' Rounds numeric values while avoiding floating-point artifacts by applying
#' a small offset before rounding.
#'
#' @param x Numeric vector.
#' @param digits Integer number of decimal places to round to.
#'
#' @return Numeric vector of rounded values.
#'
#' @details
#' A small perturbation proportional to the sign of \code{x} is added before
#' rounding to mitigate issues due to floating-point representation.
#'
#' @examples
#' round_exact(0.145, 2)
#' round_exact(c(1.005, -1.005), 2)
#'
#' @noRd
#' @keywords internal
round_exact <- function(x, digits = 0) {
  # Scale factor for the target decimal place
  scale <- 10^digits
  x_scaled <- x * scale

  # Check whether x is within floating-point noise of a half-integer
  # (i.e., the fractional part is very close to 0.5)
  frac <- abs(x_scaled) - floor(abs(x_scaled))
  near_half <- abs(frac - 0.5) < 1e-9 * abs(x_scaled + (x_scaled == 0))

  if (near_half) {
    # Nudge by a relative epsilon so the half-integer case rounds away from zero
    nudge <- sign(x) * abs(x) * 1e-10
    round(x + nudge, digits)
  } else {
    round(x, digits)
  }
}

#' Random generation from a Dirichlet distribution
#'
#' Generates random samples from a Dirichlet distribution using gamma
#' variates.
#'
#' @param n Integer number of observations to generate.
#' @param alpha Numeric vector or matrix of positive concentration parameters.
#'
#' @return A numeric matrix with \code{n} rows, where each row is a sample from
#'   the Dirichlet distribution and sums to 1.
#'
#' @details
#' Each sample is obtained by drawing independent gamma random variables and
#' normalizing them to sum to one. If \code{alpha} is a vector, it is recycled
#' across rows.
#'
#' @examples
#' rdirichlet(5, c(1, 1, 1))
#' rdirichlet(3, c(2, 5, 3))
#'
#' @export
rdirichlet <- function (n, alpha) {
  n <- as.numeric(n)
  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow = 1)
  }
  if (prod(dim(alpha)) == 0) {
    stop("alpha should be non-empty.")
  }
  if (isTRUE(any(alpha <= 0))) {
    stop("alpha must be positive.")
  }
  if (n == 1) {
    n <- nrow(alpha)
  }
  if (n > nrow(alpha)) {
    alpha <- matrix(alpha, nrow = n, ncol = ncol(alpha), byrow = TRUE)
  }
  x <- matrix(rgamma(ncol(alpha) * n, alpha), ncol = ncol(alpha))
  x/rowSums(x)
}

# -----------------------------------------------------------------------------
# Internal helper: discrete quantile lookup
#
# Given a vector of probabilities `p` and the exact CDF vector `cdf` evaluated
# over the integer support vector `supp` (sorted ascending), returns
# Q(p) = min{ x in supp : F(x) >= p } for each element of p.
# lower.tail and log.p transformations are handled by the caller before
# invoking this helper.
# -----------------------------------------------------------------------------
discrete_quantile <- function(p, cdf, supp) {
  cdf[length(cdf)] <- 1   # clamp upper endpoint for floating-point safety
  vapply(p, function(pi) {
    idx <- which(cdf >= pi)
    if (length(idx) == 0L) supp[length(supp)] else supp[idx[1L]]
  }, integer(1L))
}
