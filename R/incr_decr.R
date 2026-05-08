#' Utility function.
#'
#' This is an auxiliary function to compute the increment as the distance
#' from the case of equiprobability.
#'
#' @param m A length-one integer vector of the number of multinomial classes.
#' @param pmax A length-one numeric vector with the probability values to
#'   compute the increments for.
#' @return A numeric matrix.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{incr_2_pmax}} for the opposite calculation.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' m <- 3:50
#' pmax <- seq(0.05, 1, 0.05)
#' incr <- pmax_2_incr(m, pmax)
#' summary(as.numeric(incr))
#'
#' @export
pmax_2_incr <- function(m, pmax) {
  tmpfnc <- function(mm, pp) (mm*pp - 1)/(mm - 1)
  res <- outer(m, pmax, tmpfnc)
  res[res < 0] <- NA
  rownames(res) <- m
  colnames(res) <- pmax
  res <- drop(res)
  res
}

#' Utility function.
#'
#' This is an auxiliary function to compute the probability values from
#' an increment with respect to the case of equiprobability.
#'
#' @param m A length-one integer vector of the number of multinomial classes.
#' @param incr A length-one numeric vector with the increment values to
#'   compute the probabilities for.
#' @return A numeric matrix.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{pmax_2_incr}} for the opposite calculation.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' m <- 3:50
#' incr <- seq(0, 1, 0.01)
#' pmax <- incr_2_pmax(m, incr)
#' summary(as.numeric(pmax))
#'
#' @export
incr_2_pmax <- function(m, incr) {
  tmpfnc <- function(mm, ii) (ii + 1/(mm - 1))*(mm - 1)/mm
  res <- outer(m, incr, tmpfnc)
  rownames(res) <- m
  colnames(res) <- incr
  res <- drop(res)
  res
}

#' Utility function.
#'
#' This is an auxiliary function to compute the decrement as the distance
#' from the case of equiprobability.
#'
#' @param m A length-one integer vector of the number of multinomial classes.
#' @param pmin A length-one numeric vector with the probability values to
#'   compute the decrements for.
#' @return A numeric matrix.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{decr_2_pmin}} for the opposite calculation.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' m <- 3:50
#' pmin <- seq(0.05, 1, 0.05)
#' decr <- pmin_2_decr(m, pmin)
#' summary(as.numeric(decr))
#'
#' @export
pmin_2_decr <- function(m, pmin) {
  tmpfnc <- function(mm, pp) (1 - mm*pp)
  res <- outer(m, pmin, tmpfnc)
  res[res < 0] <- NA
  rownames(res) <- m
  colnames(res) <- pmin
  res <- drop(res)
  res
}

#' Utility function.
#'
#' This is an auxiliary function to compute the probability values from
#' an decrement with respect to the case of equiprobability.
#'
#' @param m A length-one integer vector of the number of multinomial classes.
#' @param decr A length-one numeric vector with the decrement values to
#'   compute the probabilities for.
#' @return A numeric matrix.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{pmin_2_decr}} for the opposite calculation.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' m <- 3:50
#' decr <- seq(0, 1, 0.01)
#' pmin <- decr_2_pmin(m, decr)
#' summary(as.numeric(pmin))
#'
#' @export
decr_2_pmin <- function(m, decr) {
  tmpfnc <- function(mm, dd) (1 - dd)/mm
  res <- outer(m, decr, tmpfnc)
  rownames(res) <- m
  colnames(res) <- decr
  res <- drop(res)
  res
}
