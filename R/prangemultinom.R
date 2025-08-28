#' Distribution function (CDF) of the range for a generic multinomial random
#' vector.
#'
#' This function calculates the distribution function (i.e. the probability to be
#' smaller than or equal to a given value) for the range of a generic
#' (i.e. not necessarily equiprobable) multinomial vector.
#'
#' @param x A numeric vector indicating the values of the range to compute
#'   the CDF for.
#' @param size A length-one integer specifying the total number of objects
#'   that are put into \emph{K} boxes in the typical multinomial experiment.
#' @param prob A numeric non-negative of length \emph{K} specifying the
#'   probability for the \emph{K} classes; is internally normalized to sum 1.
#' @param log A length-one logical vector; if TRUE, log probabilities are
#'   computed.
#' @param verbose A length-one logical vector; if TRUE, details are printed
#'   during the calculation.
#' @param method A length-one character vector indicating the method to use.
#'   Possible values are \code{"R"}, \code{"Rcpp"} and \code{"Corrado"}.
#' @param parallel A length-one character vector indicating the type of parallel
#'   operation to be used (if any). Possible values are \code{"multicore"}
#'   (which works only on Unix/mcOS), \code{"snow"} and \code{"no"} (i.e. serial
#'   instead of parallel computing).
#' @param threads A length-one numeric vector for the number of chains to run.
#'   If greater than 1, package \pkg{\link{parallel}} is used to take advantage of
#'   any multiprocessing or distributed computing capabilities that may be available.
#' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
#'   \code{parallel = "snow"}. If not supplied, a cluster on the local machine
#'   is created for the duration of the call.
#' @param env An optional \code{\link{environment}}. If not supplied, one is
#'   created inside the function.
#' @param tol An optional length-one non-negative value used to round the CDF
#'   values when it is close to 1.
#' @return A numeric vector providing the probabilities of the range.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{drangemulitnom}} for computing the
#'   probability mass function (PMF) of the range for a generic
#'   multinomial random vector.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   range, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' k <- 4 # dimension of the multinomial random vector
#' n <- 60
#' set.seed(101)
#' probs <- rdirichlet(1, rep(1/k, k))  # or rep(1/k, k) 
#' xseq <- 0:n
#' 
#' cdfrange <- prangemultinom(x = xseq, size = n, prob = probs, log = FALSE,
#'                            verbose = TRUE, method = "Corrado")
#' plot(xseq, cdfrange, type = "s", xlab = "x", ylab = "CDF")
#'
#' @export
prangemultinom <- function(x, size, prob, log = FALSE, verbose = FALSE, method = "Rcpp",
  parallel = "no", threads = 1, cl = NULL, env = NULL, tol = 1e-10) {
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && threads > 1L) {
    if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) {
      warning("number of cores forced to 1 (i.e. no parallel computing used).")
      threads <- 1L
    } else {
      nc <- parallel::detectCores()
      if (threads > nc) {
        warning(
          paste0("number of threads specified is larger than the number of available cores; threads set to ", nc, "."))
        threads <- nc
      }
    }
    loadNamespace("parallel") # get this out of the way before recording seed
  }
  if (!method == "Corrado") {
    stop("method not available.")
  }

  k <- length(prob)
  xlen <- length(x)
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities must be finite, non-negative and not all 0.")
  prob <- prob/s
  # x <- as.integer(x + 0.5)
  i0 <- prob == 0
  if (any(i0)) {
    prob <- prob[!i0]
    k <- length(prob)
  }
  if (is.null(env)) {
    env <- environment()
  }
  
  if (have_mc || have_snow) {
    prangemultinom_parallel <- function(x.c, size.c, prob.c, env.c, method.c) {
      if (method.c == "Corrado") {
        prangemultinom_corrado_one(x = x.c, size = size.c, prob = prob.c, verbose = FALSE, tol = tol)
      }
    }

    if (verbose) {
      devout <- ""
      if (.Platform$OS.type != "windows" && !have_mc) {
        message("--- STARTING PARALLEL CALCULATION OF ", xlen, " VALUES ---")
      } else {
        message("Performing parallel calculation of ", xlen, " values...")
      }
    } else {
      if (.Platform$OS.type != "windows") {
        devout <- '/dev/null'
      } else {
        devout <- 'nul:'
      }
    }

    res <- if (have_mc) {
             parallel::mclapply(x, function(xi) prangemultinom_parallel(x.c = xi,
               size.c = size, prob.c = prob, env.c = env, method.c = method),
               mc.cores = threads)
           } else if (have_snow) {
             if (is.null(cl)) {
               # outfile doesn't work on Windows
               cl <- parallel::makePSOCKcluster(rep("localhost", threads), outfile = devout)
             }
             parallel::setDefaultCluster(cl)
             parallel::clusterEvalQ(NULL, suppressMessages(require(XOMultinom)))
             parallel::clusterExport(NULL, c("size", "prob", "env", "method"), envir = env)
             res <- parallel::parLapply(NULL, x, function(xi) prangemultinom_parallel(x.c = xi,
               size.c = size, prob.c = prob, env.c = env, method.c = method))
             parallel::stopCluster(cl)
             res
           }
    res <- unlist(res)

    if (verbose) {
      if (.Platform$OS.type != "windows" && !have_mc){
        message("--- END OF PARALLEL CALCULATION OF ", xlen, " VALUES ---\n")
      } else {
        # message("done!")
      }
    }
  } else {
    if (method == "Corrado") {
      res <- prangemultinom_corrado(x, size, prob, log, verbose, tol)
    }
  }

  return(res)
}
