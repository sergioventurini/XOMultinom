library(XOMultinom)

### MAXIMUM ###

k <- 6 # dimension of the multinomial random vector
n <- 200
# set.seed(101)
# probs <- rdirichlet(1, rep(1/k, k))
probs <- rdirichlet(1, rep(1, k))
# probs <- rep(1/k, k)
xseq <- 0:n

# CDF
system.time(pmaxmultinom(x = 110, size = n, prob = probs, log = FALSE,
                         verbose = FALSE, method = "Rcpp"))

system.time(pmaxmultinom(x = 110, size = n, prob = probs, log = FALSE,
                         verbose = FALSE, method = "Corrado"))
