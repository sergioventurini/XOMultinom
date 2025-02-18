# Define a function to dynamically create and execute nested loops
dynamic_nested_loops_min <- function(levels, action, x, size, prob, envir = .GlobalEnv, tol) {
  # pb <- txtProgressBar(min = 0, max = (size - x + 1)^(length(prob) - 2), style = 3)
  # pb_idx <<- 1

  # Recursive function to create loops
  create_loops_min <- function(current_level, indices, xa, sz, prb, env, tol) {
    if (current_level > levels) {
      # Base case: All loops are complete, perform the action
      # setTxtProgressBar(pb, pb_idx)
      # res <- get("res", envir = env)
      # if (res > 1 - tol) {
      #   res <- 1
      #   assign("res", res, envir = env)
      #   return(invisible(NULL))
      # } else {
        action(indices, xa, env)
      # }
    } else {
      # Recursive case: Create the current loop and recurse
      for (i in xa:sz) {
        # pb_idx <<- pb_idx + 1
        create_loops_min(current_level + 1, c(indices, i), xa, sz, prb, env, tol)
      }
    }
  }
  
  # Start the recursive loop creation with level 1 and an empty index vector
  create_loops_min(1, integer(0), x, size, prob, envir, tol)
}

# Define the action to be performed inside the innermost loop
pmin_cond <- function(indices, x, envir = .GlobalEnv) {
  prx_sum <- get("prx_sum", envir = envir)
  prmin_sum <- get("prmin_sum", envir = envir)
  res <- get("res", envir = envir)

  km2 <- length(indices)
  indices_sum <- sum(indices)
  tmp <- prmin_sum[indices_sum - km2*x + 1]
  ix <- indices_sum_sub <- 0
  for (j in 1:km2) {
    prx_mat <- prx_sum[[j]]
    ix <- indices[j]
    indices_sum_sub <- ifelse(j == 1, 0, sum(indices[1:(j - 1)]))
    tmp <- tmp*prx_mat[ix - x + 1, indices_sum_sub - (j - 1)*x + 1]
  }
  res <- res + tmp

  assign("res", res, envir = envir)
}

px_cond_min <- function(x, size, prob) {
  if (size >= 0) {
    tmp <- dbinom(x, size, prob)
  } else {
    tmp <- 0
  }
  tmp
}

#' @export
pminmultinom_R_one <- function(x, size, prob, env, tol) {
  k <- length(prob)
  km2 <- k - 2
  x_tmp <- x + 1  # needed to convert P(min >= x) to P(min <= x)
  prob_c <- cumsum(prob)

  res <- 0
  assign("res", res, envir = env)

  if (x >= 0 & x < size) {
    prmin_sum <- numeric(km2*(size - x_tmp) + 1)
    for (i in (km2*x_tmp):(km2*size)) {
      if ((x_tmp <= (size - i)/2) & (x_tmp >= 0)) {
        tmp <-
          pbinom(size - i - x_tmp, size - i, prob[2]/prob_c[2]) -
          pbinom(x_tmp - 1, size - i, prob[2]/prob_c[2])
      } else if (x_tmp > (size - i)/2) {
        tmp <- 0
      } else if (x_tmp < 0) {
        tmp <- 1
      }
      prmin_sum[i - km2*x_tmp + 1] <- tmp
    }
    assign("prmin_sum", prmin_sum, envir = env)

    prx_sum <- vector(mode = "list", length = km2)
    for (p in 1:km2) {
      prx_mat <- matrix(0, nrow = (size - x_tmp + 1), ncol = ((p - 1)*(size - x_tmp) + 1))
      for (i in x_tmp:size) {
        for (j in ((p - 1)*x_tmp):((p - 1)*size)) {
          prx_mat[i - x_tmp + 1, j - (p - 1)*x_tmp + 1] <-
            px_cond_min(i, size - j, prob[k - p + 1]/prob_c[k - p + 1])
        }
      }
      prx_sum[[p]] = prx_mat
    }
    assign("prx_sum", prx_sum, envir = env)

    dynamic_nested_loops_min(km2, pmin_cond, x_tmp, size, prob, envir = env, tol)

    1 - get("res", envir = env)
  } else if (x < 0) {
    0
  } else {
    1
  }
}

#' @export
pminmultinom_R <- function(x, size, prob, log = FALSE, verbose = FALSE, env, tol) {
  xlen <- length(x)

  r <- numeric(xlen)
  for (m in 1:xlen) {
    if (verbose) print(paste0("computing P(min(X1,..., Xk) <= ", x[m], ")..."), quote = FALSE)
    if (xlen > 1 && all(diff(x) == 1) && m > 1 && abs(1 - r[m - 1]) < tol) {
      r[m] <- 1
    } else {
      r[m] <- pminmultinom_R_one(x[m], size, prob, env, tol)
    }
  }

  if (!log)
    r
  else log(r)
}

#' Distribution function (CDF) of the minimum for a generic multinomial random
#' vector.
#'
#' This function calculates the distribution function (i.e. the probability to be
#' smaller than or equal to a given value) for the minimum of a generic
#' (i.e. not necessarily equiprobable) multinomial vector.
#'
#' @param x A numeric vector indicating the values of the minimum to compute
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
#'   If greater than 1, package \pkg{\link{parallel}} is used to take advantage of any
#'   multiprocessing or distributed computing capabilities that may be available.
#' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
#'   \code{parallel = "snow"}. If not supplied, a cluster on the local machine
#'   is created for the duration of the \code{dmbc()} call.
#' @param env An optional \code{\link{environment}} . If not supplied, one is created
#'   inside the function.
#' @param tol An optional length-one non-negative value used to round the CDF
#'   values when it is close to 1.
#' @return A numeric vector providing the probabilities of the minimum.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{pminmulitnom_C}} for computing the
#'   CDF of the minimum for a generic multinomial random vector.
#' @seealso \code{\link{dminmulitnom_C}} for computing the
#'   probability mass function (PMF) of the minimum for a generic
#'   multinomial random vector.
#' @seealso \code{\link{dminmulitnom_C}} for computing the
#'   probability mass function (PMF) of the minimum for a generic
#'   multinomial random vector.
#' @references
#'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
#'   distributions of some functions of the ordered multinomial counts:
#'   minimum, minimum, range and sums of order statistics", Royal Society
#'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
#' @examples
#' k <- 4 # dimension of the multinomial random vector
#' n <- 60
#' set.seed(101)
#' probs <- rdirichlet(1, rep(1/k, k))  # or rep(1/k, k) 
#' xseq <- 0:n
#' 
#' cdfmin <- pminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
#'                        verbose = TRUE, method = "Rcpp")
#' plot(xseq, cdfmin, type = "s", xlab = "x", ylab = "CDF")
#'
#' @export
pminmultinom <- function(x, size, prob, log = FALSE, verbose = FALSE, method = "Rcpp",
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
    pminmultinom_parallel <- function(x.c, size.c, prob.c, env.c, method.c) {
      if (method.c == "Rcpp") {
        pminmultinom_C_one(x = x.c, size = size.c, prob = prob.c, verbose = FALSE, tol = tol)
      } else if (method.c == "R") {
        pminmultinom_R_one(x = x.c, size = size.c, prob = prob.c, env = env.c, tol = tol)
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
             parallel::mclapply(x, function(xi) pminmultinom_parallel(x.c = xi,
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
             res <- parallel::parLapply(NULL, x, function(xi) pminmultinom_parallel(x.c = xi,
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
      res <- pminmultinom_C(x, size, prob, log, verbose, tol)
    } else if (method == "R") {
      res <- pminmultinom_R(x, size, prob, log, verbose, env, tol)
    }
  }

  return(res)
}
