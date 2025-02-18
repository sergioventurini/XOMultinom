library(XOMultinom)

# MAX

k <- 5 # dimension of the multinomial random vector
n <- 100
# set.seed(101)
# probs <- rdirichlet(1, rep(1/k, k))
probs <- rdirichlet(1, rep(1, k))
probs <- c(0.273972229, 0.03811056154, 0.3105543744, 0.2310658149, 0.1462970201)
# probs <- rep(1/k, k)
xseq <- 0:n

# CDF
system.time(res <- pmaxmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                verbose = FALSE, method = "Rcpp", tol = 1e-7))
print(res, digits = 15)
identical(sum(diff(res)) + res[1], 1)

# MIN

k <- 5 # dimension of the multinomial random vector
n <- 100
# set.seed(101)
# probs <- rdirichlet(1, rep(1/k, k))
probs <- rdirichlet(1, rep(1, k))
probs <- c(0.273972229, 0.03811056154, 0.3105543744, 0.2310658149, 0.1462970201)
# probs <- rep(1/k, k)
xseq <- 0:n

# CDF
system.time(res <- pminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                verbose = FALSE, method = "Rcpp", tol = 1e-7))
print(res, digits = 15)
# system.time(res <- pminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
#                                 verbose = TRUE, method = "R", tol = 1e-7))
# print(res, digits = 15)
identical(sum(diff(res)) + res[1], 1)
