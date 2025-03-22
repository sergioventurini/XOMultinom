library(XOMultinom)

### MAXIMUM ###

k <- 5 # dimension of the multinomial random vector
n <- 100
# set.seed(101)
# probs <- rdirichlet(1, rep(1/k, k))
probs <- rdirichlet(1, rep(1, k))
probs <- rep(1/k, k)
xseq <- 0:n

# CDF
system.time(pmaxmultinom(x = 25, size = n, prob = probs, log = FALSE,
                         verbose = FALSE, method = "Rcpp"))

system.time(pmaxmultinom(x = 25, size = n, prob = probs, log = FALSE,
                         verbose = FALSE, method = "Corrado"))
