library(XOMultinom)

### MAXIMUM ###

k <- 4 # dimension of the multinomial random vector
n <- 100
# set.seed(101)
probs <- rdirichlet(1, rep(1/k, k))
# probs <- rep(1/k, k)
xseq <- 0:n

# CDF
cdfmax <- pmaxmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                       verbose = TRUE, method = "Rcpp",
                       parallel = "multicore", threads = 10)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq, cdfmax, type = "s", xlab = "x", ylab = "CDF", ylim = c(0, 1))

# PMF
pmfmax <- dmaxmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                       verbose = TRUE, method = "Rcpp",
                       parallel = "multicore", threads = 10)
plot(xseq, pmfmax, type = "h", xlab = "x", ylab = "PMF", ylim = c(0, 1))
points(xseq, pmfmax, pch = 20)

# dmaxmultinom(x = c(3, 10, 2, 15), size = n, prob = probs, log = FALSE,
#              verbose = TRUE, method = "Rcpp")

### Comparison with simulated data ###
set.seed(1406)
Nsim <- 1e6
simdata <- rmultinom(Nsim, n, probs)
dat_emp <- apply(simdata, 2, max)
res_emp <- ecdf(dat_emp)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(res_emp, xlim = c(0, n), main = "", xlab = "x",
     ylab = "CDF of max", ylim = c(0, 1))
lines(xseq, cdfmax, type = "s", col = "green")
dff <- res_emp(xseq) - cdfmax
summary(dff)
plot(xseq, dff, pch = 20, xlab = "x", ylab = "Difference")
lines(lowess(xseq, dff, f = 2/3), col = 2)

pmfmax_emp <- prop.table(table(dat_emp))
xidx <- match(as.integer(names(pmfmax_emp)), xseq)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq[xidx], pmfmax_emp, pch = 20, , type = "h", xlab = "x",
     ylab = "PMF of max", ylim = c(0, 1))
lines(xseq[xidx], pmfmax[xidx], pch = 20, , type = "h", col = "green")
dff <- as.numeric(pmfmax_emp) - pmfmax[xidx]
summary(dff)
plot(xseq[xidx], dff, pch = 20, xlab = "x", ylab = "PMF difference")
lines(lowess(xseq[xidx], dff, f = 2/3), col = 2)

### MINIMUM ###

library(XOMultinom)

k <- 4 # dimension of the multinomial random vector
n <- 100
# set.seed(101)
probs <- rdirichlet(1, rep(1/k, k))
# probs <- rep(1/k, k)
xseq <- 0:n

# CDF
cdfmin <- pminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                       verbose = TRUE, method = "Rcpp",
                       parallel = "multicore", threads = 10)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq, cdfmin, type = "s", xlab = "x", ylab = "CDF", ylim = c(0, 1))

# PMF
pmfmin <- dminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                       verbose = TRUE, method = "Rcpp",
                       parallel = "multicore", threads = 10)
plot(xseq, pmfmin, type = "h", xlab = "x", ylab = "PMF",
     ylim = c(0, 1))
points(xseq, pmfmin, pch = 20)

# dminmultinom(x = c(3, 10, 2, 15), size = n, prob = probs, log = FALSE,
#              verbose = TRUE, method = "Rcpp")

### Comparison with simulated data ###
set.seed(1406)
Nsim <- 1e6
simdata <- rmultinom(Nsim, n, probs)
dat_emp <- apply(simdata, 2, min)
res_emp <- ecdf(dat_emp)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(res_emp, xlim = c(0, n), main = "", xlab = "x",
     ylab = "CDF of min", ylim = c(0, 1))
lines(xseq, cdfmin, type = "s", col = "green")
dff <- res_emp(xseq) - cdfmin
summary(dff)
plot(xseq, dff, pch = 20, xlab = "x", ylab = "CDF Difference")
lines(lowess(xseq, dff, f = 2/3), col = 2)

pmfmin_emp <- prop.table(table(dat_emp))
xidx <- match(as.integer(names(pmfmin_emp)), xseq)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq[xidx], pmfmin_emp, pch = 20, , type = "h", xlab = "x",
     ylab = "PMF of min", ylim = c(0, 1))
lines(xseq[xidx], pmfmin[xidx], pch = 20, , type = "h", col = "green")
dff <- as.numeric(pmfmin_emp) - pmfmin[xidx]
summary(dff)
plot(xseq[xidx], dff, pch = 20, xlab = "x", ylab = "PMF difference")
lines(lowess(xseq[xidx], dff, f = 2/3), col = 2)
