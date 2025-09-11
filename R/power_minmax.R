power_minmax <- function(n, k, delta, power = 0.8, alpha = 0.05,
                         type, method_prob, verbose = TRUE) {
  probs_H0 <- rep(1/k, k)
  if (verbose)
    cat("  * probability % change = ", ifelse(type == "min", -delta, delta)*100, "%\n", sep = "")
  if (type == "max") {
    pmax <- incr_2_pmax(k, delta)
    probs_H1 <- c(pmax, rep((1 - pmax)/(k - 1), k - 1))
  }
  else if (type == "min") {
    pmin <- decr_2_pmin(k, delta)
    probs_H1 <- c(pmin, rep((1 - pmin)/(k - 1), k - 1))
  }
  else {
    stop("the type argument can only be equal to 'max' or 'min'.")
  }

  # k_alpha <- find_k_alpha(probs_H0, n, alpha, type = type, method_prob = method_prob)
  # gamma_prob <- find_gamma_prob(probs_H0, n, alpha, k_alpha, type = type, method_prob = method_prob)
  k_gamma <- find_k_gamma(probs_H0, n, alpha, type = type, method_prob = method_prob)
  k_alpha <- k_gamma[["k_alpha"]]
  gamma_prob <- k_gamma[["gamma_prob"]]

  if (type == "max") {
    px <- dmaxmultinom(x = (k_alpha - 1):n, size = n, prob = probs_H1,
                       log = FALSE, verbose = FALSE, method_prob = method_prob,
                       parallel = "multicore",
                       threads = parallel::detectCores(), tol = 1e-5)
    if (!is.na(gamma_prob)) {
      prb <- sum(px[-1]) + gamma_prob*px[1]
    }
  }
  else if (type == "min") {
    if (!is.na(k_alpha)) {
      Fx <- pminmultinom(x = k_alpha:(k_alpha + 1), size = n, prob = probs_H1,
                         log = FALSE, verbose = FALSE, method_prob = method_prob,
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

  if (verbose)
    cat("    - distance to power = ", power - prb, " - n value = ", n, "\n", sep = "")

  (prb - power)^2
}

power_optim <- function(n_start = 10, k, delta, power = 0.8, alpha = 0.05,
                        type, method_prob, verbose = TRUE) {
  optim(n_start, power_minmax, k = k, delta = delta, power = power, alpha = alpha,
        type = type, method_prob = method_prob, verbose = verbose)

}

power_optimize <- function(k, n_max = 5000, delta, power = 0.8, alpha = 0.05,
                           type, method_prob, verbose = TRUE) {
  optimize(power_minmax, interval = c(1, n_max), tol = 0.5,
           k = k, delta = delta, power = power, alpha = alpha, type = type,
           method_prob = method_prob, verbose = verbose)
}
