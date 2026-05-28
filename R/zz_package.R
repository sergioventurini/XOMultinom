#' XOMultinom: Exact distributions of ordered multinomial counts
#'
#' @description
#' The \pkg{XOMultinom} package provides functions for computing exact
#' distributions of selected functions of ordered multinomial counts, including
#' the maximum, minimum, range, and sums of order statistics.
#'
#' Main functions include:
#' \itemize{
#'   \item \code{\link{dmaxmultinom}()}, \code{\link{pmaxmultinom}()},
#'     \code{\link{qmaxmultinom}()}, \code{\link{rmaxmultinom}()}
#'   \item \code{\link{dminmultinom}()}, \code{\link{pminmultinom}()},
#'     \code{\link{qminmultinom}()}, \code{\link{rminmultinom}()}
#'   \item \code{\link{drangemultinom}()}, \code{\link{prangemultinom}()},
#'     \code{\link{qrangemultinom}()}, \code{\link{rrangemultinom}()}
#'   \item \code{\link{dJlargemultinom}()}, \code{\link{pJlargemultinom}()},
#'     \code{\link{qJlargemultinom}()}, \code{\link{rJlargemultinom}()}
#' }
#'
#' @name XOMultinom-package
#' @aliases XOMultinom-package XOMultinom-pkg XOMultinom
#' @docType package
#'
#' @keywords multinomial min max range
#'
#' @useDynLib XOMultinom, .registration = TRUE
#' @importFrom methods setClass setGeneric setMethod new is
#' @importFrom Rcpp evalCpp
#' @importFrom utils globalVariables packageDescription
#' @importFrom tools file_path_as_absolute
#' @importFrom graphics legend lines points segments
#' @importFrom stats dnorm median optimize pnorm rgamma uniroot
#' @importFrom rlang .data
#'
"_PACKAGE"
