library(XOMultinom)

### Section 4.3

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


### Section 5.1.5

library(ggplot2)
library(patchwork)

### Maximum test ###

pow_seq <- c(0.8, 0.9)
alpha_seq <- c(0.05, 0.01, 0.001)
m_seq <- c(3:10, 15, 20, 30)
incr_seq <- seq(0.1, 0.9, 0.1)
n_master <- list()
for (pow in pow_seq) {
  for (alpha in alpha_seq) {
    cat("Power = ", pow, " - alpha = ", alpha, "\n", sep = "")
    res <- maxmin_multinom_size(m_seq, incr_seq, power = pow, alpha = alpha,
                                n_max = 200, type = "max",
                                verbose = TRUE, optmethod = "uniroot")
    n_master[[paste0("power = ", pow, " - alpha = ", alpha)]] <- res$sizes
  }
}

df <- data.frame(power = numeric(0),
                 alpha = numeric(0),
                 m = numeric(0),
                 n = numeric(0),
                 incr = numeric(0))
i <- 1
for (p in 1:length(pow_seq)) {
  for (a in 1:length(alpha_seq)) {
    for (m in 1:length(m_seq)) {
      df_tmp <- data.frame(power = pow_seq[p],
                           alpha = alpha_seq[a],
                           m = m_seq[m],
                           n = n_master[[i]][[m]],
                           incr = incr_seq*100)
      df <- rbind(df, df_tmp)
    }
    i <- i + 1
  }
}
df$alpha <- factor(df$alpha,
                   levels = sort(unique(df$alpha), decreasing = TRUE))

p_max <- ggplot(df, aes(x = incr, y = n, color = factor(m), group = m)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.5) +
  facet_grid(alpha ~ power, labeller = label_both, scales = "free_y") +
  scale_color_viridis_d(option = "turbo", name = "m") +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  labs(
    x = "% increase in one probability",
    y = "Sample size"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1))
p_max <- p_max +
  labs(title = "Maximum test") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  
### Minimum test ###

pow_seq <- c(0.8, 0.9)
alpha_seq <- c(0.05, 0.01, 0.001)
m_seq <- c(3:10)
decr_seq <- seq(0.2, 0.9, 0.1)
n_master <- list()
for (pow in pow_seq) {
  for (alpha in alpha_seq) {
    cat("Power = ", pow, " - alpha = ", alpha, "\n", sep = "")
    res <- maxmin_multinom_size(m_seq, decr_seq, power = pow, alpha = alpha,
                                n_max = 1000, type = "min",
                                verbose = TRUE, optmethod = "uniroot")
    n_master[[paste0("power = ", pow, " - alpha = ", alpha)]] <- res$sizes
  }
}

df <- data.frame(power = numeric(0),
                 alpha = numeric(0),
                 m = numeric(0),
                 n = numeric(0),
                 decr = numeric(0))
i <- 1
for (p in 1:length(pow_seq)) {
  for (a in 1:length(alpha_seq)) {
    for (m in 1:length(m_seq)) {
      df_tmp <- data.frame(power = pow_seq[p],
                           alpha = alpha_seq[a],
                           m = m_seq[m],
                           n = n_master[[i]][[m]],
                           decr = decr_seq*100)
      df <- rbind(df, df_tmp)
    }
    i <- i + 1
  }
}
df$alpha <- factor(df$alpha,
                   levels = sort(unique(df$alpha), decreasing = TRUE))

p_min <- ggplot(df, aes(x = decr, y = n, color = factor(m), group = m)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.5) +
  facet_grid(alpha ~ power, labeller = label_both, scales = "free_y") +
  scale_color_viridis_d(option = "turbo", name = "m") +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  labs(
    x = "% decrease in one probability",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 1))
p_min <- p_min +
  labs(title = "Minimum test") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

p_max + p_min
