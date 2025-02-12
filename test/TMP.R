library(XOMultinom)

# MAX

k <- 3 # dimension of the multinomial random vector
n <- 5
set.seed(101)
probs <- rdirichlet(1, rep(1/k, k))
probs <- rep(1/k, k)
xseq <- 0:n

# CDF
system.time(res <- pmaxmultinom(x = 4, size = n, prob = probs, log = FALSE,
                    verbose = TRUE, method = "Rcpp", tol = 1e-7))
res

# MIN

k <- 5 # dimension of the multinomial random vector
n <- 20
set.seed(101)
probs <- rdirichlet(1, rep(1/k, k))
probs <- rep(1/k, k)
xseq <- 0:n

# CDF
system.time(res <- pminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                verbose = FALSE, method = "Rcpp", tol = 1e-7))
res
system.time(res <- pminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                verbose = FALSE, method = "R", tol = 1e-7))
res
