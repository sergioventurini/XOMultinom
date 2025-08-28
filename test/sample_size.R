library(XOMultinom)

# MAXIMUM

pow_seq <- c(0.8, 0.9)
alpha_seq <- c(0.05, 0.01, 0.001)
k_seq <- c(3:10, 15, 20, 30)
incr_seq <- seq(0.1, 0.9, 0.1)
n_seq <- numeric(length(incr_seq))
n_master <- list()
for (pow in pow_seq) {
  for (alpha in alpha_seq) {
    cat("Power = ", pow, " - alpha = ", alpha, "\n", sep = "")
    n_all <- maxmin_multinom_size(k_seq, incr_seq, power = pow, alpha = alpha,
                                  type = "max", method = "Corrado", verbose = FALSE)
    n_master[[paste0("power = ", pow, " - alpha = ", alpha)]] <- n_all
  }
}

save(n_master, file = "/Users/Sergio/dev/XOMultinom/data/sample_size_MAX.RData")

# toplot <- 3
# plot(incr_seq, n_master[[toplot]][[1]], type = "n", xlim = range(incr_seq),
#      ylim = c(0, max(unlist(n_master[[toplot]]))))
# for (k in 1:length(n_master[[toplot]])) {
#   lines(incr_seq, n_master[[toplot]][[k]], type = "b", lwd = 2, col = k)
# }

###

# MINIMUM

pow_seq <- c(0.8, 0.9)
alpha_seq <- c(0.05, 0.01, 0.001)
k_seq <- c(3:10, 15, 20, 30)
decr_seq <- seq(0.1, 0.9, 0.1)
n_seq <- numeric(length(decr_seq))
n_master <- list()
for (pow in pow_seq) {
  for (alpha in alpha_seq) {
    cat("Power = ", pow, " - alpha = ", alpha, "\n", sep = "")
    n_all <- maxmin_multinom_size(k_seq, decr_seq, power = pow, alpha = alpha,
                                  type = "min", method = "Corrado", verbose = FALSE)
    n_master[[paste0("power = ", pow, " - alpha = ", alpha)]] <- n_all
  }
}

toplot <- 3
plot(incr_seq, n_master[[toplot]][[1]], type = "n", xlim = range(incr_seq),
     ylim = c(0, max(unlist(n_master[[toplot]]))))
for (k in 1:length(n_master[[toplot]])) {
  lines(incr_seq, n_master[[toplot]][[k]], type = "b", lwd = 2, col = k)
}

save(n_master, file = "/Users/Sergio/dev/XOMultinom/data/sample_size_MIN.RData")

### TMP ###

# pow_seq <- 0.8
# alpha_seq <-0.05
# k_seq <- 40
# incr_seq <- 0.9
# probs_H0 <- rep(1/k_seq, k_seq)

# debugonce(find_k_gamma)
# nval <- 4
# find_k_gamma(probs_H0, nval, alpha_seq, type = "max", method = "Corrado")

# debugonce(maxmin_multinom_size)
# maxmin_multinom_size(k_seq, incr_seq, power = pow_seq, alpha = alpha_seq,
#                      type = "max", method = "Corrado", verbose = TRUE)
