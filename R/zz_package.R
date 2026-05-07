#' XOMultinom: Exact distributions of ordered multinomial counts
#'
#' @description
#' The \pkg{XOMultinom} package provides functions for computing exact
#' distributions of selected functions of ordered multinomial counts, including
#' the maximum, minimum, range, and sums of order statistics.
#'
#' Main functions include:
#' \itemize{
#'   \item \code{\link{dmaxmultinom}()}
#'   \item \code{\link{dminmultinom}()}
#'   \item \code{\link{drangemultinom}()}
#'   \item \code{\link{pmaxmultinom}()}
#'   \item \code{\link{pminmultinom}()}
#'   \item \code{\link{prangemultinom}()}
#' }
#'
#' @name XOMultinom-package
#' @aliases XOMultinom-package XOMultinom-pkg XOMultinom
#' @docType package
#'
#' @useDynLib XOMultinom, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils globalVariables packageDescription
#' @importFrom tools file_path_as_absolute
#'
"_PACKAGE"
