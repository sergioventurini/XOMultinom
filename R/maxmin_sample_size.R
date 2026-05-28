#' Critical value and randomization probability for max/min tests
#'
#' Computes the critical value \eqn{k_\alpha} and the corresponding
#' randomization probability \eqn{\gamma} for hypothesis tests based on the
#' maximum or minimum of a multinomial random vector under the null hypothesis
#' of equiprobable categories.
#'
#' @param probs Numeric vector of probabilities. Must correspond to the
#'   equiprobable case (i.e., all equal).
#' @param n Integer number of trials.
#' @param alpha Significance level in (0, 1).
#' @param type Character string; either \code{"max"} or \code{"min"} indicating
#'   the test statistic.
#'
#' @return A list with components:
#'   \item{k_alpha}{Critical value.}
#'   \item{gamma_prob}{Randomization probability.}
#'
#' @details
#' The function determines the rejection region for tests based on the
#' maximum or minimum cell count. When the test is not exact, a randomized
#' decision rule is constructed via \eqn{\gamma}.
#'
#' @export
find_k_gamma <- function(probs, n, alpha = 0.05, type) {
  if (alpha <= 0 | alpha >= 1)
    stop("the alpha argument must be in between 0 and 1.")
  if (any(probs != probs[1]))
    stop("the probabilities must correspond to the equiprobability case.")

  if (type == "max") {

    xval <- 0:(n + 1)
    px   <- dmaxmultinom(x = xval, size = n, prob = probs,
                         log = FALSE, verbose = FALSE)$values

    pgtk    <- rev(cumsum(rev(px)))
    k_alpha <- min(xval[which(pgtk <= alpha)])

    idx           <- match(k_alpha - 1, xval)
    p_k_alpha     <- sum(px[(idx + 1):length(px)])
    p_k_alpha_m_1 <- px[idx]

    gamma_prob <- if (p_k_alpha_m_1 > 0)
                    (alpha - p_k_alpha) / p_k_alpha_m_1
                  else
                    0

  } else if (type == "min") {

    # binary search
    k <- length(probs)

    # Step 1 - does the test have any rejection region at all?
    # P(min <= 0 | H0) should be > 0 for any finite n, but check defensively.
    F0 <- pminmultinom(x = 0L, size = n, prob = probs,
                       log = FALSE, verbose = FALSE)$values
    if (is.na(F0) || F0 > alpha) {
      return(list(k_alpha = NA, gamma_prob = NA))
    }

    # Step 2 - find an upper bound hi such that F(hi) > alpha
    # Start at floor(n/k) (the mean of each cell) and double if needed
    hi <- max(1L, floor(n / k))
    Fhi <- pminmultinom(x = hi, size = n, prob = probs,
                        log = FALSE, verbose = FALSE)$values
    while (!is.na(Fhi) && Fhi <= alpha && hi < n) {
      hi  <- min(n - 1L, hi * 2L)
      Fhi <- pminmultinom(x = hi, size = n, prob = probs,
                          log = FALSE, verbose = FALSE)$values
    }
    if (is.na(Fhi) || Fhi <= alpha) {
      return(list(k_alpha = NA, gamma_prob = NA))
    }

    # Step 3 - binary search in [lo=0, hi]
    lo <- 0L
    while (hi - lo > 1L) {
      mid  <- (lo + hi) %/% 2L
      Fmid <- pminmultinom(x = mid, size = n, prob = probs,
                           log = FALSE, verbose = FALSE)$values
      if (is.na(Fmid)) {
        # Computation failed mid-search; use the last safe lower bound
        return(list(k_alpha = lo, gamma_prob = 0))
      }
      if (Fmid <= alpha) lo <- mid else hi <- mid
    }
    k_alpha <- lo

    if (k_alpha >= n - 1L) {
      return(list(k_alpha = k_alpha, gamma_prob = 0))
    }

    # Step 4 - compute gamma from the two flanking CDF values
    Flo <- pminmultinom(x = k_alpha,       size = n, prob = probs,
                        log = FALSE, verbose = FALSE)$values
    Fup <- pminmultinom(x = k_alpha + 1L,  size = n, prob = probs,
                        log = FALSE, verbose = FALSE)$values

    if (is.na(Flo) || is.na(Fup)) {
      return(list(k_alpha = k_alpha, gamma_prob = 0))
    }

    p_mass <- Fup - Flo
    gamma_prob <- if (p_mass > .Machine$double.eps)
                    (alpha - Flo) / p_mass
                  else
                    0
    gamma_prob <- min(max(gamma_prob, 0), 1)  # clamp fp noise

  } else {
    stop("the type argument can only be equal to 'max' or 'min'.")
  }

  list(k_alpha = k_alpha, gamma_prob = gamma_prob)
}

#' Critical value for max/min multinomial tests
#'
#' Computes the critical value \eqn{k_\alpha} for hypothesis tests based on
#' the maximum or minimum of a multinomial random vector.
#'
#' @param probs Numeric vector of probabilities. Must correspond to the
#'   equiprobable case.
#' @param n Integer number of trials.
#' @param alpha Significance level in (0, 1).
#' @param type Character string; either \code{"max"} or \code{"min"}.
#'
#' @return Integer critical value \eqn{k_\alpha}. Returns \code{NA} if no valid
#'   rejection region exists.
#'
#' @export
find_k_alpha <- function(probs, n, alpha = 0.05, type) {
  if (alpha <= 0 | alpha >= 1)
    stop("the alpha argument must be in between 0 and 1.")
  if (any(probs != probs[1]))
    stop("the probabilities must correspond to the equiprobability case.")

  if (type == "max") {
    xval <- 0:(n + 1)
    px <- dmaxmultinom(x = xval, size = n, prob = probs,
                       log = FALSE, verbose = FALSE)$values
    pgtk <- rev(cumsum(rev(px)))
    k_alpha <- min(xval[which(pgtk <= alpha)])
  }
  else if (type == "min") {
    xval <- 0:n
    Fx <- pminmultinom(x = xval, size = n, prob = probs,
                       log = FALSE, verbose = FALSE)$values
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

#' Randomization probability for max/min multinomial tests
#'
#' Computes the randomization probability \eqn{\gamma} associated with a
#' critical value \eqn{k_\alpha} for tests based on the maximum or minimum
#' of a multinomial random vector.
#'
#' @param probs Numeric vector of probabilities. Must correspond to the
#'   equiprobable case.
#' @param n Integer number of trials.
#' @param alpha Significance level in (0, 1).
#' @param k_alpha Integer critical value.
#' @param type Character string; either \code{"max"} or \code{"min"}.
#'
#' @return Numeric value representing the randomization probability. Returns
#'   \code{NA} if not defined.
#'
#' @export
find_gamma_prob <- function(probs, n, alpha = 0.05, k_alpha, type) {
  if (alpha <= 0 | alpha >= 1)
    stop("the alpha argument must be in between 0 and 1.")
  if (any(probs != probs[1]))
    stop("the probabilities must correspond to the equiprobability case.")

  if (type == "max") {
    px <- dmaxmultinom(x = (k_alpha - 1):n, size = n, prob = probs,
                       log = FALSE, verbose = FALSE)$values
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
                         log = FALSE, verbose = FALSE)$values
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

#' Sample size determination for multinomial max/min tests
#'
#' Computes the required sample size to achieve a target power for hypothesis
#' tests based on the maximum or minimum of a multinomial random vector under
#' deviations from equiprobability.
#'
#' @param m_seq Integer vector of numbers of categories.
#' @param change_seq Numeric vector of probability perturbations from the
#'   equiprobable case.
#' @param power Desired power level in (0, 1).
#' @param alpha Significance level in (0, 1).
#' @param n_max Maximum sample size considered in the search.
#' @param type Character string; either \code{"max"} or \code{"min"}.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed.
#' @param optmethod Character string; optimization method, either
#'   \code{"uniroot"} or \code{"optimize"}.
#' @param extendInt Passed to \code{uniroot()} when used.
#'
#' @return A list where each element corresponds to a value of \code{m_seq}
#'   and contains the required sample sizes for each value in \code{change_seq}.
#'
#' @details
#' The function evaluates the sample size needed to detect deviations from
#' equiprobability with a given power, using tests based on either the
#' maximum or minimum multinomial cell count.
#'
#' @examples
#' \donttest{
#' pow <- 0.8
#' alpha <- 0.05
#' m_seq <- 3:8
#' incr_seq <- seq(0.2, 0.8, 0.1)
#' res <- maxmin_multinom_size(m_seq, incr_seq, power = pow, alpha = alpha,
#'                             n_max = 200, type = "max",
#'                             verbose = TRUE, optmethod = "uniroot")
#' summary(res)
#' plot(res)
#' }
#'
#' @export
maxmin_multinom_size <- function(m_seq, change_seq, power = 0.8, alpha = 0.05,
                                 n_max = 500, type, verbose = TRUE,
                                 optmethod = "uniroot", extendInt = "upX") {
  if (power <= 0 | power >= 1)
    stop("the power argument must be in between 0 and 1.")
  if (alpha <= 0 | alpha >= 1)
    stop("the alpha argument must be in between 0 and 1.")
  if (type != "max" && type != "min")
    stop("the type argument can only be equal to 'max' or 'min'.")

  n_seq <- numeric(length(change_seq))
  n_all <- list()
  for (m in m_seq) {
    if (verbose)
      cat("  Number of classes = ", m, "\n", sep = "")

    i <- 1
    for (change in change_seq) {
      if (verbose)
        cat("  * probability % change = ", ifelse(type == "min", -change, change)*100, "%\n", sep = "")

      if (optmethod == "optimize") {
        res <- power_optimize(m, n_max = n_max, change = change, power = power, alpha = alpha,
                              type = type, verbose = verbose)
        n_max_new <- n_max
        if (type == "max") {
          while (res$objective > 1e-3) {
            if (verbose)
              message("  minimization algorithm stuck --> restarting reducing the maximum n value")
            n_max_new <- ceiling(n_max_new*.9)
            if (n_max_new < 2)
              stop("the minimization algorithm did not converge --> try increasing the maximum n value")
            res <- power_optimize(m, n_max = n_max_new, change = change, power = power, alpha = alpha,
                                  type = type, verbose = verbose)
          }
        }
        else if (type == "min") {
          while (res$objective > 1e-3) {
            if (verbose)
              message("  minimization algorithm stuck --> restarting increasing the maximum n value")
            n_max_new <- ceiling(n_max_new*1.1)
            if (n_max_new > 2500)
              stop("the minimization algorithm did not converge --> try increasing the maximum n value")
            res <- power_optimize(m, n_max = n_max_new, change = change, power = power, alpha = alpha,
                                  type = type, verbose = verbose)
          }
        }

        n_seq[i] <- ceiling(res$minimum)
      }
      else if (optmethod == "uniroot") {
        res <- power_uniroot(m, n_max = n_max, change = change, power = power, alpha = alpha,
                            type = type, extendInt = extendInt, verbose = verbose)
        
        n_seq[i] <- ceiling(res$root)
      }

      i <- i + 1
    }
    names(n_seq) <- change_seq
    n_all[[paste0("m = ", m)]] <- n_seq
  }

  return(new_xomultinom_size(sizes      = n_all,
                             m_seq      = m_seq,
                             change_seq = change_seq,
                             power      = power,
                             alpha      = alpha,
                             type       = type))
}
