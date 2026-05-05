library(XOMultinom)

### MAXIMUM ###

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

# CDF
system.time(cdfmax <- pmaxmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                   verbose = TRUE))
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq, cdfmax$values, type = "s", xlab = "x", ylab = "CDF", ylim = c(0, 1))

# PMF
system.time(pmfmax <- dmaxmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                   verbose = TRUE))
plot(xseq, pmfmax$values, type = "h", xlab = "x", ylab = "PMF")
points(xseq, pmfmax$values, pch = 20)

### Comparison with simulated data ###
set.seed(1406)
Nsim <- 1e6
simdata <- rmultinom(Nsim, n, probs)
dat_emp <- apply(simdata, 2, max)
res_emp <- ecdf(dat_emp)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(res_emp, xlim = c(0, n), main = "", xlab = "x",
     ylab = "CDF of max", ylim = c(0, 1), pch = 20)
lines(xseq, cdfmax$values, type = "s", col = "green")
dff <- res_emp(xseq) - cdfmax$values
summary(dff)
plot(xseq, dff, pch = 20, xlab = "x", ylab = "Difference")
lines(lowess(xseq, dff, f = 2/3), col = 2)

pmfmax_emp <- prop.table(table(dat_emp))
xidx <- match(as.integer(names(pmfmax_emp)), xseq)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq[xidx], pmfmax_emp, pch = 20, , type = "h", xlab = "x",
     ylab = "PMF of max")
lines(xseq[xidx], pmfmax$values[xidx], pch = 20, , type = "h", col = "green")
dff <- as.numeric(pmfmax_emp) - pmfmax$values[xidx]
summary(dff)
plot(xseq[xidx], dff, pch = 20, xlab = "x", ylab = "PMF difference")
lines(lowess(xseq[xidx], dff, f = 2/3), col = 2)

### MINIMUM ###

library(XOMultinom)

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

# CDF
system.time(cdfmin <- pminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                   verbose = TRUE))
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq, cdfmin$values, type = "s", xlab = "x", ylab = "CDF", ylim = c(0, 1))

# PMF
system.time(pmfmin <- dminmultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                   verbose = TRUE))
plot(xseq, pmfmin$values, type = "h", xlab = "x", ylab = "PMF")
points(xseq, pmfmin$values, pch = 20)

### Comparison with simulated data ###
set.seed(1406)
Nsim <- 1e6
simdata <- rmultinom(Nsim, n, probs)
dat_emp <- apply(simdata, 2, min)
res_emp <- ecdf(dat_emp)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(res_emp, xlim = c(0, n), main = "", xlab = "x",
     ylab = "CDF of min", ylim = c(0, 1), pch = 20)
lines(xseq, cdfmin$values, type = "s", col = "green")
dff <- res_emp(xseq) - cdfmin$values
summary(dff)
plot(xseq, dff, pch = 20, xlab = "x", ylab = "CDF Difference")
lines(lowess(xseq, dff, f = 2/3), col = 2)

pmfmin_emp <- prop.table(table(dat_emp))
xidx <- match(as.integer(names(pmfmin_emp)), xseq)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq[xidx], pmfmin_emp, pch = 20, , type = "h", xlab = "x",
     ylab = "PMF of min")
lines(xseq[xidx], pmfmin$values[xidx], pch = 20, , type = "h", col = "green")
dff <- as.numeric(pmfmin_emp) - pmfmin$values[xidx]
summary(dff)
plot(xseq[xidx], dff, pch = 20, xlab = "x", ylab = "PMF difference")
lines(lowess(xseq[xidx], dff, f = 2/3), col = 2)

### RANGE ###

library(XOMultinom)

k <- 8 # dimension of the multinomial random vector
n <- 50
# set.seed(101)
# probs <- rdirichlet(1, rep(1/k, k))
probs <- rdirichlet(1, rep(1, k))
# probs <- rep(1/k, k)
xseq <- 0:n

# CDF
system.time(cdfrange <- prangemultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                       verbose = TRUE))
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq, cdfrange$values, type = "s", xlab = "x", ylab = "CDF", ylim = c(0, 1))

# PMF
system.time(pmfrange <- drangemultinom(x = xseq, size = n, prob = probs, log = FALSE,
                                       verbose = TRUE))
plot(xseq, pmfrange$values, type = "h", xlab = "x", ylab = "PMF")
points(xseq, pmfrange$values, pch = 20)

### Comparison with simulated data ###
set.seed(1406)
Nsim <- 1e6
simdata <- rmultinom(Nsim, n, probs)
dat_emp <- apply(simdata, 2, function(x) diff(range(x)))
res_emp <- ecdf(dat_emp)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(res_emp, xlim = c(0, n), main = "", xlab = "x",
     ylab = "CDF of range", ylim = c(0, 1), pch = 20)
lines(xseq, cdfrange$values, type = "s", col = "green")
dff <- res_emp(xseq) - cdfrange$values
summary(dff)
plot(xseq, dff, pch = 20, xlab = "x", ylab = "CDF Difference")
lines(lowess(xseq, dff, f = 2/3), col = 2)

pmfrange_emp <- prop.table(table(dat_emp))
xidx <- match(as.integer(names(pmfrange_emp)), xseq)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1) + 0.1)
plot(xseq[xidx], pmfrange_emp, pch = 20, , type = "h", xlab = "x",
     ylab = "PMF of range")
lines(xseq[xidx], pmfrange$values[xidx], pch = 20, , type = "h", col = "green")
dff <- as.numeric(pmfrange_emp) - pmfrange$values[xidx]
summary(dff)
plot(xseq[xidx], dff, pch = 20, xlab = "x", ylab = "PMF difference")
lines(lowess(xseq[xidx], dff, f = 2/3), col = 2)
