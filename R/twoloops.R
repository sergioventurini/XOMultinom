twoloops_R <- function(d, n) {
  s <- numeric(d)
  np <- n^(0:d)
  k <- np[d + 1]
  
  for (u in 0:(k - 1)) {
    for (v in 0:(d - 1)) {
      tmp <- 1 + (u %% np[v + 2]) %/% np[v + 1]
      s[v + 1] <- tmp
      # cat(s[v + 1], " ")
    }
    # cat("\n")
  }
}

twoloops_matrix_R <- function(d, n) {
  np <- n^(0:d)
  k <- np[d + 1]
  
  u_vals <- 0:(k - 1)
  indices <- outer(u_vals, np[2:(d + 1)], function(x, y) 1 + (x %% y) %/% (y / n))
  
  # print(tail(indices))
}
