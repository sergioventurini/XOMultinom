# Define a function to dynamically create and execute nested loops
dynamic_nested_loops_max <- function(levels, action, x, size, prob, envir = .GlobalEnv) {
  # Recursive function to create loops
  create_loops <- function(current_level, indices, xa, sz, prb, env) {
    if (current_level > levels) {
      # Base case: All loops are complete, perform the action
      action(indices, env)
    } else {
      # Recursive case: Create the current loop and recurse
      for (i in 0:xa) {
        create_loops(current_level + 1, c(indices, i), xa, sz, prb, env)
      }
    }
  }
  
  # Start the recursive loop creation with level 1 and an empty index vector
  create_loops(1, integer(0), x, size, prob, envir)
}

# Define the action to be performed inside the innermost loop
pmax_cond <- function(indices, envir = .GlobalEnv) {
  prx_sum <- get("prx_sum", envir = envir)
  prmax_sum <- get("prmax_sum", envir = envir)
  res <- get("res", envir = envir)

  km2 <- length(indices)
  indices_sum <- sum(indices)
  tmp <- prmax_sum[indices_sum + 1]
  ix <- indices_sum_sub <- 0
  for (j in 1:km2) {
    prx_mat <- prx_sum[[j]]
    ix <- indices[j]
    indices_sum_sub <- ifelse(j == 1, 0, sum(indices[1:(j - 1)]))
    tmp <- tmp*prx_mat[ix + 1, indices_sum_sub + 1]
  }
  res <- res + tmp

  assign("res", res, envir = envir)
}

px_cond_max <- function(x, size, prob) {
  if (size >= 0) {
    tmp <- dbinom(x, size, prob)
  } else {
    tmp <- 0
  }
  tmp
}

#' @export
pmaxmultinom_R_one <- function(x, size, prob, env) {
  k <- length(prob)
  km2 <- k - 2
  prob_c <- cumsum(prob)

  res <- 0
  assign("res", res, envir = env)

  if (x >= 0 & x < size) {
    prmax_sum <- numeric(km2*x + 1)
    for (i in 0:(km2*x)) {
      if ((x <= (size - i)) & (x >= (size - i)/2)) {
        tmp <- 
          pbinom(x, size - i, prob[2]/prob_c[2]) -
          pbinom(size - i - x - 1, size - i, prob[2]/prob_c[2])
      }
      else if (x > (size - i)) {
        tmp <- 1
      }
      else if (x < (size - i)/2) {
        tmp <- 0
      }
      prmax_sum[i + 1] <- tmp
    }
    assign("prmax_sum", prmax_sum, envir = env)

    prx_sum <- vector(mode = "list", length = km2)
    for (p in 1:km2) {
      prx_mat <- matrix(0, nrow = (x + 1), ncol = ((p - 1)*x + 1))
      for (i in 0:x) {
        for (j in 1:((p - 1)*x + 1)) {
          prx_mat[i + 1, j] <- px_cond_max(i, size - j + 1, prob[k - p + 1]/prob_c[k - p + 1])
        }
      }
      prx_sum[[p]] = prx_mat
    }
    assign("prx_sum", prx_sum, envir = env)

    dynamic_nested_loops_max(km2, pmax_cond, x, size, prob, envir = env)

    get("res", envir = env)
  } else if (x < 0) {
    0
  } else {
    1
  }
}

#' @export
pmaxmultinom_R <- function(x, size, prob, log = FALSE, verbose = FALSE, env, tol) {
  xlen <- length(x)

  r <- numeric(xlen)
  for (m in 1:xlen) {
    if (verbose) print(paste0("computing P(max(X1,..., Xk) <= x) for x = ", x[m], "..."), quote = FALSE)
    if (xlen > 1 && all(diff(x) == 1) && m > 1 && abs(1 - r[m - 1]) < tol) {
      r[m] <- 1
    } else {
      r[m] <- pmaxmultinom_R_one(x[m], size, prob, env)
    }
  }

  if (!log) 
    r
  else log(r)
}

#' Distribution function (CDF) of the maximum for a generic multinomial random
#' vector.
#'
#' This function calculates the distribution function (i.e. the probability to be
#' smaller than or equal to a given value) for the maximum of a generic
#' (i.e. not necessarily equiprobable) multinomial vector.
#'
#' @param x A numeric vector indicating the values of the maximum to compute
#'   the CDF for.
#' @param size A length-one integer specifying the total number of objects
#'   that are put into \emph{K} boxes in the typical multinomial experiment.
#' @param prob A numeric non-negative of length \emph{K} specifying the
#'   probability for the \emph{K} classes; is internally normalized to sum 1.
#' @param logd A length-one logical vector; if TRUE, log probabilities are
#'   computed.
#' @param verbose A length-one logical vector; if TRUE, details are printed
#'   during the calculation.
#' @param parallel A length-one character vector indicating the type of parallel
#'   operation to be used (if any). Possible values are \code{multicore}
#'   (which worksonly on Unix/mcOS), \code{snow} and \code{no} (i.e. serial
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
#' @return A numeric vector providing the probabilities of the maximum.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{pminmulitnom_C}} for computing the
#'   CDF of the minimum for a generic multinomial random vector.
#' @seealso \code{\link{dmaxmulitnom_C}} for computing the
#'   probability mass function (PMF) of the maximum for a generic
#'   multinomial random vector.
#' @seealso \code{\link{dminmulitnom_C}} for computing the
#'   probability mass function (PMF) of the minimum for a generic
#'   multinomial random vector.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   maximum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' k <- 4 # dimension of the multinomial random vector
#' n <- 60
#' set.seed(101)
#' probs <- rdirichlet(1, rep(1/k, k))  # or rep(1/k, k) 
#' xseq <- 0:n
#' 
#' cdfmax <- pmaxmultinom(x = xseq, size = n, prob = probs, log = FALSE,
#'                        verbose = TRUE, method = "Rcpp")
#' plot(xseq, cdfmax, type = "s", xlab = "x", ylab = "CDF")
#'
#' @export
pmaxmultinom <- function(x, size, prob, log = FALSE, verbose = FALSE, method = "Rcpp",
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
  if (!method == "R" && !method == "Rcpp") {
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
    pmaxmultinom_parallel <- function(x.c, size.c, prob.c, env.c, method.c) {
      if (method.c == "Rcpp") {
        pmaxmultinom_C_one(x = x.c, size = size.c, prob = prob.c, this_env = env.c)
      } else if (method.c == "R") {
        pmaxmultinom_R_one(x = x.c, size = size.c, prob = prob.c, env = env.c)
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
             parallel::mclapply(x, function(xi) pmaxmultinom_parallel(x.c = xi,
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
             res <- parallel::parLapply(NULL, x, function(xi) pmaxmultinom_parallel(x.c = xi,
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
    if (method == "Rcpp") {
      res <- pmaxmultinom_C(x, size, prob, log, verbose, env, tol)
    } else if (method == "R") {
      res <- pmaxmultinom_R(x, size, prob, log, verbose, env, tol)
    }
  }

  return(res)
}
