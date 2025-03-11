library(XOMultinom)

m <- 25 # dimension of the multinomial random vector
n <- 100

probs_null <- rep(1/m, m)
r <- 11
alpha <- 1 - prangemultinom(x = r - 1, size = n, prob = probs_null, log = FALSE,
                            verbose = FALSE, method = "Corrado")

set.seed(1406)
nsim <- 5e5

simdata <- rmultinom(nsim, n, probs_null)
chisq_vals <- apply(simdata, 2,
                    function(x) chisq.test(as.table(x), p = probs_null)$statistic)
sum(chisq_vals > qchisq(p = 1 - alpha, df = m - 1) - m/n)/nsim

###

pi_prime <- .20
probs <- c(pi_prime, rep((1 - pi_prime)/(m - 1), m - 1))
r <- 11
1 - prangemultinom(x = r - 1, size = n, prob = probs, log = FALSE,
                   verbose = FALSE, method = "Corrado")
simdata <- rmultinom(nsim, n, probs)
chisq_vals <- apply(simdata, 2,
                    function(x) chisq.test(as.table(x), p = probs_null)$statistic)
sum(chisq_vals > qchisq(p = 1 - alpha, df = m - 1) - m/n)/nsim

###

# section 5.3
m <- 9 # dimension of the multinomial random vector
n <- 25
probs_null <- rep(1/m, m)
pminmultinom(x = 0, size = n, prob = probs_null, log = FALSE,
               verbose = FALSE, method = "Corrado")
