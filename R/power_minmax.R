power_minmax <- function(n, m, change, power = 0.8, alpha = 0.05,
                         type, verbose = TRUE, optmethod) {
  probs_H0 <- rep(1/m, m)

  if (type == "max") {
    pmax <- incr_2_pmax(m, change)
    probs_H1 <- c(pmax, rep((1 - pmax)/(m - 1), m - 1))
  }
  else if (type == "min") {
    pmin <- decr_2_pmin(m, change)
    probs_H1 <- c(pmin, rep((1 - pmin)/(m - 1), m - 1))
  }

  # k_alpha <- find_k_alpha(probs_H0, n, alpha, type = type)
  # gamma_prob <- find_gamma_prob(probs_H0, n, alpha, k_alpha, type = type)
  k_gamma <- find_k_gamma(probs_H0, n, alpha, type = type)
  k_alpha <- k_gamma[["k_alpha"]]
  gamma_prob <- k_gamma[["gamma_prob"]]

  prb <- 0
  if (type == "max") {
    px <- dmaxmultinom(x = (k_alpha - 1):n, size = n, prob = probs_H1,
                       log = FALSE, verbose = FALSE)$values
    if (!is.na(gamma_prob)) {
      prb <- sum(px[-1]) + gamma_prob*px[1]
    }
  }
  else if (type == "min") {
    if (!is.na(k_alpha)) {
      Fx <- pminmultinom(x = k_alpha:(k_alpha + 1), size = n, prob = probs_H1,
                         log = FALSE, verbose = FALSE)$values
      if (!is.na(gamma_prob)) {
        prb <- Fx[1] + gamma_prob*diff(Fx)
      }
    }
  }

  if (verbose)
    cat("    - distance to power = ", power - prb, " - n value = ", n, "\n", sep = "")

  if (optmethod == "optimize") {
    (prb - power)^2
  }
  else if (optmethod == "uniroot") {
    prb - power
  }
}

power_optimize <- function(m, n_max = 3000, change, power = 0.8, alpha = 0.05,
                           type, verbose = TRUE) {
  optimize(power_minmax, interval = c(1, n_max), maximum = FALSE, tol = 0.5,
           m = m, change = change, power = power, alpha = alpha, type = type,
           verbose = verbose, optmethod = "optimize")
}

power_uniroot <- function(m, n_max = 500, change, power = 0.8, alpha = 0.05,
                           type, extendInt = "upX",
                           verbose = TRUE) {
  uniroot(power_minmax, lower = 1, upper = n_max, extendInt = extendInt,
          tol = 0.5, maxiter = 10000,
          m = m, change = change, power = power, alpha = alpha, type = type,
          verbose = verbose, optmethod = "uniroot")
}
