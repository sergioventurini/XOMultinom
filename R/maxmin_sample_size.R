#' Utility function.
#'
#' This is an auxiliary function to compute the increment as the distance
#' from the case of equiprobability.
#'
#' @param k A length-one integer vector of the number of multinomial classes.
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
#' k <- 3:50
#' pmax <- seq(0.05, 1, 0.05)
#' incr <- pmax_2_incr(k, pmax)
#' summary(as.numeric(incr))
#'
#' @export
find_k_gamma <- function(probs, n, alpha = 0.05, type, method) {
  if (alpha <= 0 | alpha >= 1)
    stop("the alpha argument must be in between 0 and 1.")
  if (any(diff(probs) > 0))
    stop("the probabilities must correspond to the equiprobability case.")

  if (type == "max") {
    xval <- 0:(n + 1)
    px <- dmaxmultinom(x = xval, size = n, prob = probs,
                       log = FALSE, verbose = FALSE, method = method,
                       parallel = "multicore",
                       threads = parallel::detectCores(), tol = 1e-5)

    pgtk <- rev(cumsum(rev(px)))
    k_alpha <- min(xval[which(pgtk <= alpha)])

    idx <- match(k_alpha - 1, xval)
    p_k_alpha <- sum(px[(idx + 1):length(px)])
    p_k_alpha_m_1 <- px[idx]
    if (p_k_alpha_m_1 > 0) {
      gamma_prob <- (alpha - p_k_alpha)/p_k_alpha_m_1
    }
    else {
      gamma_prob <- NA
    }
  }
  else if (type == "min") {
    xval <- 0:n
    Fx <- pminmultinom(x = xval, size = n, prob = probs,
                       log = FALSE, verbose = FALSE, method = method,
                       parallel = "multicore",
                       threads = parallel::detectCores(), tol = 1e-5)
    chk <- which(Fx <= alpha)
    if (length(chk) > 0) {
      k_alpha <- max(xval[chk])

      idx <- match(k_alpha, xval)
      p_k_alpha <- Fx[idx]
      p_k_alpha_p_1 <- diff(Fx[c(idx, idx + 1)])
      if (p_k_alpha_p_1 > 0) {
        gamma_prob <- (alpha - p_k_alpha)/p_k_alpha_p_1
      }
      else {
        gamma_prob <- NA
      }
    }
    else {
      k_alpha <- NA
      gamma_prob <- NA
    }

  }
  else {
    stop("the type argument can only be equal to 'max' or 'min'.")
  }

  list(k_alpha = k_alpha, gamma_prob = gamma_prob)
}

#' Utility function.
#'
#' This is an auxiliary function to compute the increment as the distance
#' from the case of equiprobability.
#'
#' @param k A length-one integer vector of the number of multinomial classes.
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
#' k <- 3:50
#' pmax <- seq(0.05, 1, 0.05)
#' incr <- pmax_2_incr(k, pmax)
#' summary(as.numeric(incr))
#'
#' @export
find_k_alpha <- function(probs, n, alpha = 0.05, type, method) {
  if (alpha <= 0 | alpha >= 1)
    stop("the alpha argument must be in between 0 and 1.")
  if (any(diff(probs) > 0))
    stop("the probabilities must correspond to the equiprobability case.")

  if (type == "max") {
    xval <- 0:(n + 1)
    px <- dmaxmultinom(x = xval, size = n, prob = probs,
                       log = FALSE, verbose = FALSE, method = method,
                       parallel = "multicore",
                       threads = parallel::detectCores(), tol = 1e-5)
    pgtk <- rev(cumsum(rev(px)))
    k_alpha <- min(xval[which(pgtk <= alpha)])
  }
  else if (type == "min") {
    xval <- 0:n
    Fx <- pminmultinom(x = xval, size = n, prob = probs,
                       log = FALSE, verbose = FALSE, method = method,
                       parallel = "multicore",
                       threads = parallel::detectCores(), tol = 1e-5)
    chk <- which(Fx <= alpha)
    if (length(chk) > 0) {
      k_alpha <- max(xval[chk])
    }
    else {
      k_alpha <- NA
    }
  }
  else {
    stop("the type argument can only be equal to 'max' or 'min'.")
  }

  k_alpha
}

#' Utility function.
#'
#' This is an auxiliary function to compute the increment as the distance
#' from the case of equiprobability.
#'
#' @param k A length-one integer vector of the number of multinomial classes.
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
#' k <- 3:50
#' pmax <- seq(0.05, 1, 0.05)
#' incr <- pmax_2_incr(k, pmax)
#' summary(as.numeric(incr))
#'
#' @export
find_gamma_prob <- function(probs, n, alpha = 0.05, k_alpha, type, method) {
  if (alpha <= 0 | alpha >= 1)
    stop("the alpha argument must be in between 0 and 1.")
  if (any(diff(probs) > 0))
    stop("the probabilities must correspond to the equiprobability case.")

  if (type == "max") {
    px <- dmaxmultinom(x = (k_alpha - 1):n, size = n, prob = probs,
                       log = FALSE, verbose = FALSE, method = method,
                       parallel = "multicore",
                       threads = parallel::detectCores(), tol = 1e-5)
    p_k_alpha <- sum(px[-1])
    p_k_alpha_m_1 <- px[1]
    if (p_k_alpha_m_1 > 0) {
      gamma_prob <- (alpha - p_k_alpha)/p_k_alpha_m_1
    }
    else {
      gamma_prob <- NA
    }
  }
  else if (type == "min") {
    if (is.na(k_alpha)) {
      gamma_prob <- NA
    }
    else {
      Fx <- pminmultinom(x = k_alpha:(k_alpha + 1), size = n, prob = probs,
                         log = FALSE, verbose = FALSE, method = method,
                         parallel = "multicore",
                         threads = parallel::detectCores(), tol = 1e-5)
      p_k_alpha <- Fx[1]
      p_k_alpha_p_1 <- diff(Fx)
      if (p_k_alpha_p_1 > 0) {
        gamma_prob <- (alpha - p_k_alpha)/p_k_alpha_p_1
      }
      else {
        gamma_prob <- NA
      }
    }
  }
  else {
    stop("the type argument can only be equal to 'max' or 'min'.")
  }

  gamma_prob
}

#' Utility function.
#'
#' This is an auxiliary function to compute the increment as the distance
#' from the case of equiprobability.
#'
#' @param k A length-one integer vector of the number of multinomial classes.
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
#' k <- 3:50
#' pmax <- seq(0.05, 1, 0.05)
#' incr <- pmax_2_incr(k, pmax)
#' summary(as.numeric(incr))
#'
#' @export
maxmin_multinom_size_OLD <- function(k_seq, change_seq, power = 0.8, alpha = 0.05,
                                     type, method, verbose = TRUE) {
  n_seq <- numeric(length(change_seq))
  n_all <- list()
  for (k in k_seq) {
    if (verbose)
      cat("  Number of classes = ", k, "\n", sep = "")
    probs_H0 <- rep(1/k, k)
    i <- 1
    for (change in change_seq) {
      if (verbose)
        cat("  * probability % change = ", change*100, "%\n", sep = "")
      if (type == "max") {
        pmax <- incr_2_pmax(k, change)
        probs_H1 <- c(pmax, rep((1 - pmax)/(k - 1), k - 1))
      }
      else if (type == "min") {
        pmin <- decr_2_pmin(k, change)
        probs_H1 <- c(pmin, rep((1 - pmin)/(k - 1), k - 1))
      }
      else {
        stop("the type argument can only be equal to 'max' or 'min'.")
      }
      
      nval <- 1
      prb <- 0
      while (prb < power) {
        if (verbose)
          cat("    - distance to power = ", power - prb, " - n value = ", nval, "\n", sep = "")
        # k_alpha <- find_k_alpha(probs_H0, nval, alpha, type = type, method = method)
        # gamma_prob <- find_gamma_prob(probs_H0, nval, alpha, k_alpha, type = type, method = method)
        k_gamma <- find_k_gamma(probs_H0, nval, alpha, type = type, method = method)
        k_alpha <- k_gamma[["k_alpha"]]
        gamma_prob <- k_gamma[["gamma_prob"]]

        if (type == "max") {
          px <- dmaxmultinom(x = (k_alpha - 1):nval, size = nval, prob = probs_H1,
                             log = FALSE, verbose = FALSE, method = method,
                             parallel = "multicore",
                             threads = parallel::detectCores(), tol = 1e-5)
          if (!is.na(gamma_prob)) {
            prb <- sum(px[-1]) + gamma_prob*px[1]
          }
        }
        else if (type == "min") {
          if (!is.na(k_alpha)) {
            Fx <- pminmultinom(x = k_alpha:(k_alpha + 1), size = nval, prob = probs_H1,
                               log = FALSE, verbose = FALSE, method = method,
                               parallel = "multicore",
                               threads = parallel::detectCores(), tol = 1e-5)
            if (!is.na(gamma_prob)) {
              prb <- Fx[1] + gamma_prob*diff(Fx)
            }
          }
        }
        else {
          stop("the type argument can only be equal to 'max' or 'min'.")
        }
        nval <- nval + 1
      }
      n_seq[i] <- nval - 1
      i <- i + 1
    }
    names(n_seq) <- change_seq
    n_all[[paste0("k = ", k)]] <- n_seq
  }

  n_all
}

#' Utility function.
#'
#' This is an auxiliary function to compute the increment as the distance
#' from the case of equiprobability.
#'
#' @param k A length-one integer vector of the number of multinomial classes.
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
#' k <- 3:50
#' pmax <- seq(0.05, 1, 0.05)
#' incr <- pmax_2_incr(k, pmax)
#' summary(as.numeric(incr))
#'
#' @export
maxmin_multinom_size <- function(k_seq, change_seq, power = 0.8, alpha = 0.05,
                                 n_max = 500, type, method, verbose = TRUE,
                                 optmethod = "uniroot", extendInt = "upX") {
  if (power <= 0 | power >= 1)
    stop("the power argument must be in between 0 and 1.")
  if (alpha <= 0 | alpha >= 1)
    stop("the alpha argument must be in between 0 and 1.")
  if (type != "max" && type != "min")
    stop("the type argument can only be equal to 'max' or 'min'.")

  n_seq <- numeric(length(change_seq))
  n_all <- list()
  for (k in k_seq) {
    if (verbose)
      cat("  Number of classes = ", k, "\n", sep = "")

    i <- 1
    for (change in change_seq) {
      if (verbose)
        cat("  * probability % change = ", ifelse(type == "min", -change, change)*100, "%\n", sep = "")

      if (optmethod == "optimize") {
        res <- power_optimize(k, n_max = n_max, change = change, power = power, alpha = alpha,
                              type = type, method_f = method, verbose = verbose)
        n_max_new <- n_max
        if (type == "max") {
          while (res$objective > 1e-3) {
            if (verbose)
              message("  minimization algorithm stuck --> restarting reducing the maximum n value")
            n_max_new <- ceiling(n_max_new*.9)
            if (n_max_new < 2)
              stop("the minimization algorithm did not converge --> try increasing the maximum n value")
            res <- power_optimize(k, n_max = n_max_new, change = change, power = power, alpha = alpha,
                                  type = type, method_f = method, verbose = verbose)
          }
        }
        else if (type == "min") {
          while (res$objective > 1e-3) {
            if (verbose)
              message("  minimization algorithm stuck --> restarting increasing the maximum n value")
            n_max_new <- ceiling(n_max_new*1.1)
            if (n_max_new > 2500)
              stop("the minimization algorithm did not converge --> try increasing the maximum n value")
            res <- power_optimize(k, n_max = n_max_new, change = change, power = power, alpha = alpha,
                                  type = type, method_f = method, verbose = verbose)
          }
        }

        n_seq[i] <- ceiling(res$minimum)
      }
      else if (optmethod == "uniroot") {
        res <- power_uniroot(k, n_max = n_max, change = change, power = power, alpha = alpha,
                            type = type, method_f = method, extendInt = extendInt, verbose = verbose)
        
        n_seq[i] <- ceiling(res$root)
      }

      i <- i + 1
    }
    names(n_seq) <- change_seq
    n_all[[paste0("k = ", k)]] <- n_seq
  }

  n_all
}
