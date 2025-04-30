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

###

k <- 10
nvals <- c(100, 150, 200, 300, 500)
times <- numeric(length(nvals))
set.seed(101)
# probs <- rdirichlet(1, rep(1/k, k))
probs <- rdirichlet(1, rep(1, k))
# probs <- rep(1/k, k)
for (i in seq_along(nvals)) {
  cat(paste0("n value = ", nvals[i], "\n"))
  xvals <- 0:nvals[i]
  times[i] <- system.time(pmaxmultinom(x = xvals, size = nvals[i], prob = probs, log = FALSE,
                                       verbose = FALSE, method = "Corrado"))[3]
}
times
