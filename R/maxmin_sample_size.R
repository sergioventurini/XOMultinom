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
maxmin_multinom_size <- function(k_seq, delta_seq, power = 0.8, alpha = 0.05,
                                 type, method, verbose = TRUE) {
  n_seq <- numeric(length(delta_seq))
  n_all <- list()
  for (k in k_seq) {
    if (verbose)
      cat("  Number of classes = ", k, "\n", sep = "")

    i <- 1
    for (delta in delta_seq) {
      res <- power_optimize(k, n_max = 5000, delta, power = power, alpha = alpha,
                            type = type, method_prob = method, verbose = verbose)

      n_seq[i] <- ceiling(res$minimum)
      i <- i + 1
    }
    names(n_seq) <- delta_seq
    n_all[[paste0("k = ", k)]] <- n_seq
  }

  n_all
}
