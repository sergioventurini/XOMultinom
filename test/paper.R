library(XOMultinom)

m <- 6
n <- 200
prob_equi <- rep(1/m, m)           # equiprobable: Bonetti et al. back-end
prob_gen  <- c(0.3, 0.15, 0.15,
               0.15, 0.15, 0.1)    # non-equiprobable: Corrado back-end

# PMF and CDF (automatic dispatch)
pmf_max <- dmaxmultinom(x = 0:n, size = n, prob = prob_equi)
cdf_max <- pmaxmultinom(x = 0:n, size = n, prob = prob_gen)

# Quantiles: median and 95th percentile
qmaxmultinom(p = c(0.5, 0.95), size = n, prob = prob_equi)

# Built-in plot and summary via S3 methods
par(mfrow = c(1, 2))
plot(pmf_max, main = "Multinomial Max - PMF - Equiprobable",
     xlab = expression(N["<1>"]), xlim = c(30, 60), cex.main = 0.85)
plot(cdf_max, main = "Multinomial Max - CDF - Non-equiprobable",
     xlab = expression(N["<1>"]), xlim = c(30, 100), cex.main = 0.85)

summary(pmf_max)

# 100,000 random draws
set.seed(1406)
n_sim <- 100000
sims <- rmaxmultinom(n = n_sim, size = n, prob = prob_equi)
par(mfrow = c(1, 1))
hist(
  sims,
  col = "#4C78A8",
  border = "white",
  lwd = 1.2,
  freq = FALSE,
  xlab = expression(N["<1>"]),
  ylab = "Density",
  main = "Random Draws for the Multinomial Maximum - Equiprobable Case",
  xlim = c(30, 55),
  ylim = c(0, 0.15),
  xaxt = "n",
  yaxt = "n"
)
axis(1, at = seq(30, 55, 5), pos = 0)
axis(2, at = seq(0, 0.15, 0.025), pos = 30, las = 1)
