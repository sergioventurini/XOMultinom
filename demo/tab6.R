### NOTE: the following code will take several minutes to run! ###

library(XOMultinom)

# Table 6 - Distribution of the multinomial minimum.
t_min <- 0
t_max <- 4 # the original t_max value was 7
t <- t_min:t_max
m <- c(5, 10, 15, 20)
n_5 <- c(seq(5, 30, 5), 40)
n_10 <- c(seq(30, 60, 10), 80)
n_15 <- 90
n_20 <- 100
n_list <- list(n_5, n_10, n_15, n_20)
names(n_list) <- m
n_mat <- matrix(NA, nrow = max(sapply(n_list, length)), ncol = length(m))
for (j in 1:ncol(n_mat)) {
  n_mat[1:length(n_list[[j]]), j] <- n_list[[j]]
}
# colnames(n_mat) <- paste0("m = ", m)
# rownames(n_mat) <- rep("n = ", nrow(n_mat))

probs <- array(NA, dim = c(nrow(n_mat), length(t), ncol(n_mat)))
dimnames(probs)[[2]] <- t
dimnames(probs)[[3]] <- m
for (m_idx in 1:length(m)) {
  for (n_idx in 1:length(n_list[[m_idx]])) {
    cat(paste0("--> m = ", m[m_idx], " - n = ", n_list[[m_idx]][n_idx], "\n"))
    cat("    t = ")
    for (t_idx in 1:length(t)) {
      if (t_idx < length(t)) {
        cat(paste0(t[t_idx], "..."))
      } else {
        cat(t[t_idx])
      }
      # probs[n_idx, t_idx, m_idx] <- 1 - smallest_order_value(t[t_idx] + 1, n_list[[m_idx]][n_idx], m[m_idx])
      probs[n_idx, t_idx, m_idx] <- 1 - smallest_order_value_C(t[t_idx] + 1, n_list[[m_idx]][n_idx], m[m_idx])
    }
    cat("\n")
  }
}

round_exact(probs, digits = 3)
