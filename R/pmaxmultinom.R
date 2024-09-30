# Define a function to dynamically create and execute nested loops
dynamic_nested_loops_max <- function(levels, action, x, size, prob, envir = .GlobalEnv) {
  # Recursive function to create loops
  create_loops <- function(current_level, indices, xa, sz, prb, env) {
    if (current_level > levels) {
      # Base case: All loops are complete, perform the action
      action(indices, xa, sz, prb, env)
    } else {
      # Recursive case: Create the current loop and recurse
      prb_c <- cumsum(prb)
      k <- length(prb)
      idx <- k - current_level + 1
      for (i in 0:xa) {
      	px <- px_cond_max(i, sz - sum(indices), prb[idx]/prb_c[idx])
        prx <- get("prx", envir = env)
        prx[[current_level]] <- c(prx[[current_level]], px)
        assign("prx", prx, envir = env)
        create_loops(current_level + 1, c(indices, i), xa, sz, prb, env)
      }
    }
  }
  
  # Start the recursive loop creation with level 1 and an empty index vector
  create_loops(1, integer(0), x, size, prob, envir)
}

# Define the action to be performed inside the innermost loop
pmax_cond <- function(indices, x, size, prob, envir = .GlobalEnv) {
  prob_c <- cumsum(prob)
  if ((x <= (size - sum(indices))) & (x >= (size - sum(indices))/2)) {
    tmp <- 
      pbinom(x, size - sum(indices), prob[2]/prob_c[2]) -
      pbinom(size - sum(indices) - x - 1, size - sum(indices), prob[2]/prob_c[2])
  }
  else if (x > (size - sum(indices))) {
    tmp <- 1
  }
  else if (x < (size - sum(indices))/2) {
    tmp <- 0
  }
  prmax <- get("prmax", envir = envir)
  # prmax <- rbind(prmax, c(indices, tmp))  # old implementation that took too much time
  prmax_idx <- get("prmax_idx", envir = envir)
  prmax[prmax_idx, ] <- c(indices, tmp)
  prmax_idx <- prmax_idx + 1
  assign("prmax", prmax, envir = envir)
  assign("prmax_idx", prmax_idx, envir = envir)
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
  if (x >= 0 & x < size) {
    prmax_idx <- 1
    assign("prmax_idx", prmax_idx, envir = env)
    prmax <- data.frame(matrix(NA, nrow = (x + 1)^(k - 2), ncol = (k - 1)))
    assign("prmax", prmax, envir = env)
    # prmax <- data.frame()  # old implementation that took too much time
    prx <- vector(mode = "list", length = k - 2)
    assign("prx", prx, envir = env)
    dynamic_nested_loops_max(k - 2, pmax_cond, x, size, prob, envir = env)

    prmax <- get("prmax", envir = env)
    prx <- get("prx", envir = env)
    red <- prmax
    idx <- k - 2
    while (idx > 0) {
      red[, ncol(red)] <- red[, ncol(red)] * prx[[idx]]
      if (idx > 1) {
        which_cols <- 1:(ncol(red) - 2)
        by_cols <- red[which_cols]
      } else {
        by_cols <- list(rep(0, nrow(red)))
      }
      red <- aggregate(red[, ncol(red)], by_cols, sum)
      red[1:(ncol(red) - 1)] <- red[rev(1:(ncol(red) - 1))]
      idx <- idx - 1
    }

    red[1, 2]
  } else if (x < 0) {
    0
  } else {
    1
  }
}

#' @export
pmaxmultinom_R <- function(x, size, prob, log = FALSE, verbose = FALSE, env) {
  xlen <- length(x)

  r <- numeric(xlen)
  for (m in 1:xlen) {
    if (verbose) print(paste0("computing P(max(X1,..., Xk) <= x) for x = ", x[m], "..."), quote = FALSE)
    r[m] <- pmaxmultinom_R_one(x[m], size, prob, env)
    gc()
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
#'   If greater than 1, package \pkg{\link{parallel}} is used to take advantage of any
#'   multiprocessing or distributed computing capabilities that may be available.
#' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
#'   \code{parallel = "snow"}. If not supplied, a cluster on the local machine
#'   is created for the duration of the call.
#' @param env An optional \code{\link{environment}} . If not supplied, one is created
#'   inside the function.
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
  parallel = "no", threads = 1, cl = NULL, env = NULL) {
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
      res <- pmaxmultinom_C(x, size, prob, log, verbose, env)
    } else if (method == "R") {
      res <- pmaxmultinom_R(x, size, prob, log, verbose, env)
    }
  }

  return(res)
}
