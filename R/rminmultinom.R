#' Random generation from the distribution of the multinomial minimum
#'
#' Draws independent random samples from the exact distribution of the minimum
#' cell count of a multinomial random vector with arbitrary cell probabilities.
#'
#' @param n Integer number of random samples to draw.
#' @param size Integer number of trials in each multinomial experiment.
#' @param prob Numeric vector of non-negative cell probabilities. Values are
#'   internally normalized to sum to 1. Categories with zero probability are
#'   removed before computation.
#'
#' @return Integer vector of length \code{n} containing independent draws from
#'   the distribution of \eqn{\min(N_1, \ldots, N_m)}.
#'
#' @details
#'   The exact PMF over the support \eqn{\{0, 1, \ldots, \mathtt{size}\}} is
#'   computed once using \code{\link{dminmultinom}}, and \code{n} independent
#'   draws are then obtained via \code{\link[base]{sample}} with those
#'   probabilities as weights.  The cost is therefore dominated by the single
#'   PMF evaluation and is independent of \code{n}.
#'
#' @examples
#' m <- 4; n <- 60
#' probs <- rep(1 / m, m)
#'
#' set.seed(42)
#' sims <- rminmultinom(n = 1000, size = n, prob = probs)
#' hist(sims, breaks = 20, main = "Simulated multinomial minima")
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' Corrado, C.J. (2011). The exact distribution of the maximum, minimum and
#' the range of Multinomial/Dirichlet and Multivariate Hypergeometric
#' frequencies. Statistical Computing, 21, 349--359.
#' \doi{10.1007/s11222-010-9174-3}
#'
#' @seealso
#' \code{\link{dminmultinom}} for the PMF,
#' \code{\link{pminmultinom}} for the CDF,
#' \code{\link{qminmultinom}} for quantiles.
#'
#' @export
rminmultinom <- function(n, size, prob) {
  supp <- 0L:size
  pmf  <- dminmultinom(x = supp, size = size, prob = prob,
                       log = FALSE, verbose = FALSE)$values
  pmf  <- pmax(pmf, 0)
  pmf  <- pmf / sum(pmf)
  sample(supp, size = n, replace = TRUE, prob = pmf)
}
