### WARNING: the following code may take several minutes to run! ###

library(XOMultinom)

# Table 5 - Distribution of the sum of the first three largest order statistics.
t_min <- 8
t_max <- 15
t <- c(t_min:t_max, seq(16, 34, by = 2))
m <- c(10, 15, 20, 25, 30, 40, 50)
n_10 <- seq(10, 60, 10)
n_15 <- seq(15, 90, 15)
n_20 <- seq(20, 100, 20)
n_25 <- seq(25, 50, 25)
n_30 <- seq(30, 60, 30)
n_40 <- seq(40, 80, 40)
n_50 <- seq(50, 100, 50)
n_list <- list(n_10, n_15, n_20, n_25, n_30, n_40, n_50)
names(n_list) <- m
n_mat <- matrix(NA, nrow = max(sapply(n_list, length)), ncol = length(m))
for (j in 1:ncol(n_mat)) {
  n_mat[1:length(n_list[[j]]), j] <- n_list[[j]]
}
# colnames(n_mat) <- paste0("m = ", m)
# rownames(n_mat) <- rep("n = ", nrow(n_mat))
J <- 3

probs <- array(NA, dim = c(nrow(n_mat), length(t), ncol(n_mat)))
dimnames(probs)[[2]] <- t
dimnames(probs)[[3]] <- m
for (m_idx in 1:length(m)) {
  for (n_idx in 1:length(n_list[[m_idx]])) {
    cat(paste0("--> m = ", m[m_idx], " - n = ", n_list[[m_idx]][n_idx], "\n"))
    cellp <- rep(1/m[m_idx], m[m_idx])
    probs[n_idx, , m_idx] <- pJlargemultinom(t, n_list[[m_idx]][n_idx],
                                             cellp, J = 3, FALSE, FALSE)$values
  }
}

round_exact(probs, digits = 3)
