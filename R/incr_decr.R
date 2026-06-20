#' Utility function.
#'
#' This is an auxiliary function to compute the increment as the distance
#' from the case of equiprobability.
#'
#' @param m An integer vector of numbers of multinomial classes.
#' @param pmax A numeric vector of probability values for which to compute
#'   the increments.
#' @return A numeric matrix.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{incr_2_pmax}} for the opposite calculation.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, \doi{10.1098/rsos.190198}.
#' @examples
#' m <- 3:50
#' pmax <- seq(0.05, 1, 0.05)
#' incr <- pmax_2_incr(m, pmax)
#' summary(as.numeric(incr))
#'
#' @noRd
#' @keywords internal
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
#' @param m An integer vector of numbers of multinomial classes.
#' @param incr A numeric vector of increment values for which to compute
#'   the probabilities.
#' @return A numeric matrix.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{pmax_2_incr}} for the opposite calculation.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, \doi{10.1098/rsos.190198}.
#' @examples
#' m <- 3:50
#' incr <- seq(0, 1, 0.01)
#' pmax <- incr_2_pmax(m, incr)
#' summary(as.numeric(pmax))
#'
#' @noRd
#' @keywords internal
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
#' @param m An integer vector of numbers of multinomial classes.
#' @param pmin A numeric vector of probability values for which to compute
#'   the decrements.
#' @return A numeric matrix.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{decr_2_pmin}} for the opposite calculation.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, \doi{10.1098/rsos.190198}.
#' @examples
#' m <- 3:50
#' pmin <- seq(0.05, 1, 0.05)
#' decr <- pmin_2_decr(m, pmin)
#' summary(as.numeric(decr))
#'
#' @noRd
#' @keywords internal
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
#' @param m An integer vector of numbers of multinomial classes.
#' @param decr A numeric vector of decrement values for which to compute
#'   the probabilities.
#' @return A numeric matrix.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{pmin_2_decr}} for the opposite calculation.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, \doi{10.1098/rsos.190198}.
#' @examples
#' m <- 3:50
#' decr <- seq(0, 1, 0.01)
#' pmin <- decr_2_pmin(m, decr)
#' summary(as.numeric(pmin))
#'
#' @noRd
#' @keywords internal
decr_2_pmin <- function(m, decr) {
  tmpfnc <- function(mm, dd) (1 - dd)/mm
  res <- outer(m, decr, tmpfnc)
  rownames(res) <- m
  colnames(res) <- decr
  res <- drop(res)
  res
}
