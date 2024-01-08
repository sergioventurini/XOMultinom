cat("NOTE: the following code will take several minutes to run!\n")

library(XOMultinom)

# Table 3 - Distribution of the multinomial maximum.
t <- 2:16
m <- c(5, 10, 15, 20, 25, 30, 40, 50)
n_5 <- seq(5, 30, 5)
n_10 <- seq(10, 60, 10)
n_15 <- seq(15, 90, 15)
n_20 <- seq(20, 120, 20)
n_25 <- seq(25, 150, 25)
n_30 <- seq(30, 180, 30)
n_40 <- seq(40, 200, 40)
n_50 <- seq(50, 200, 50)
n_list <- list(n_5, n_10, n_15, n_20, n_25, n_30, n_40, n_50)
names(n_list) <- m
n_mat <- matrix(NA, nrow = max(sapply(n_list, length)), ncol = length(m))
for (j in 1:ncol(n_mat)) {
  n_mat[1:length(n_list[[j]]), j] <- n_list[[j]]
}
# colnames(n_mat) <- paste0("m = ", m)
# rownames(n_mat) <- rep("n = ", nrow(n_mat))

probs <- array(NA, dim = c(length(t), dim(n_mat)))
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
      probs[t_idx, n_idx, m_idx] <- max_order_statistic(t[t_idx], n_list[[m_idx]][n_idx], m[m_idx])
    }
    cat("\n")
  }
}

