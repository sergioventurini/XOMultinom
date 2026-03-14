library(microbenchmark)

k <- 10 # dimension of the multinomial random vector
n <- 300
# set.seed(101)
# probs <- rdirichlet(1, rep(1/k, k))
probs <- rdirichlet(1, rep(1, k))
# probs <- rep(1/k, k)
r <- 0.97
probs <- r^(0:(k - 1))
probs <- probs/sum(probs)
xseq <- 0:n

microbenchmark(
  max = pmaxmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                     verbose = FALSE),
  min = pminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                     verbose = FALSE),
  times = 10
)
