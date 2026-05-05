#' XOMultinom: Exact distributions of ordered multinomial counts
#'
#' @description
#' The \pkg{XOMultinom} package provides functions for computing exact
#' distributions of selected functions of ordered multinomial counts, including
#' the maximum, minimum, range, and sums of order statistics.
#'
#' Main functions include:
#' \itemize{
#'   \item \code{\link{highest_order_statistics}()}
#'   \item \code{\link{max_order_statistic}()}
#'   \item \code{\link{range_probability}()}
#'   \item \code{\link{smallest_order_value}()}
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
#' @importFrom ggplot2 autoplot
#' @importFrom Rcpp evalCpp
#' @importFrom utils globalVariables packageDescription
#' @importFrom tools file_path_as_absolute
#'
"_PACKAGE"
