#' Random generation from the distribution of the sum of \eqn{J} largest
#' order statistics for a multinomial distribution
#'
#' Draws independent random samples from the exact distribution of
#' \eqn{S_J = \sum_{j=1}^{J} N_{\langle j \rangle}} for a multinomial random
#' vector with equal cell probabilities.
#'
#' @param n Integer number of random samples to draw.
#' @param size Integer number of trials in each multinomial experiment.
#' @param prob Numeric vector of non-negative, equal cell probabilities.
#'   Only the equiprobable case is supported; a non-equiprobable \code{prob}
#'   will raise an error (propagated from \code{\link{dJlargemultinom}}).
#' @param J Integer number of largest order statistics to sum. Defaults to
#'   \code{2}.
#'
#' @return Integer vector of length \code{n} containing independent draws from
#'   the distribution of \eqn{S_J}.
#'
#' @details
#'   The exact PMF over the support \eqn{\{0, 1, \ldots, \mathtt{size}\}} is
#'   computed once using \code{\link{dJlargemultinom}}, and \code{n} independent
#'   draws are then obtained via \code{\link[base]{sample}} with those
#'   probabilities as weights.  Only the equiprobable case is supported,
#'   consistent with \code{\link{dJlargemultinom}}.
#'
#' @examples
#' m <- 4; n <- 60
#' probs <- rep(1 / m, m)
#'
#' set.seed(42)
#' sims <- rJlargemultinom(n = 1000, size = n, prob = probs, J = 3)
#' hist(sims, breaks = 20, main = "Simulated sums of 3 largest order statistics")
#'
#' @references
#' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
#' distributions of some functions of the ordered multinomial counts:
#' maximum, minimum, range and sums of order statistics. Royal Society
#' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
#'
#' @seealso
#' \code{\link{dJlargemultinom}} for the PMF,
#' \code{\link{pJlargemultinom}} for the CDF,
#' \code{\link{qJlargemultinom}} for quantiles.
#'
#' @export
rJlargemultinom <- function(n, size, prob, J = 2) {
  supp <- 0L:size
  pmf  <- dJlargemultinom(x = supp, size = size, prob = prob, J = J,
                          log = FALSE, verbose = FALSE)$values
  pmf  <- pmax(pmf, 0)
  pmf  <- pmf / sum(pmf)
  sample(supp, size = n, replace = TRUE, prob = pmf)
}
