library(XOMultinom)
library(ggplot2)

# MAXIMUM

pow_seq <- c(0.8, 0.9)
alpha_seq <- c(0.05, 0.01, 0.001)
k_seq <- c(3:10, 15, 20, 30)
incr_seq <- seq(0.1, 0.9, 0.1)
n_seq <- numeric(length(incr_seq))
n_master <- list()
for (pow in pow_seq) {
  for (alpha in alpha_seq) {
    cat("Power = ", pow, " - alpha = ", alpha, "\n", sep = "")
    n_all <- maxmin_multinom_size(k_seq, incr_seq, power = pow, alpha = alpha,
                                  type = "max", method = "Corrado", verbose = FALSE)
    n_master[[paste0("power = ", pow, " - alpha = ", alpha)]] <- n_all
  }
}

save(n_master, file = "/Users/Sergio/dev/XOMultinom/data/sample_size_MAX.RData")

pow_seq <- c(0.8, 0.9)
alpha_seq <- c(0.05, 0.01, 0.001)
k_seq <- c(3:10, 15, 20, 30)
incr_seq <- seq(0.1, 0.9, 0.1)
df <- data.frame(power = numeric(0),
                 alpha = numeric(0),
                 k = numeric(0),
                 n = numeric(0),
                 incr = numeric(0))
i <- 1
for (p in 1:length(pow_seq)) {
  for (a in 1:length(alpha_seq)) {
    for (k in 1:length(k_seq)) {
      df_tmp <- data.frame(power = pow_seq[p],
                           alpha = alpha_seq[a],
                           k = k_seq[k],
                           n = n_master[[i]][[k]],
                           incr = incr_seq*100)
      df <- rbind(df, df_tmp)
    }
    i <- i + 1
  }
}
df$alpha <- factor(df$alpha, 
                   levels = sort(unique(df$alpha), decreasing = TRUE))

ggplot(df, aes(x = incr, y = n, color = factor(k), group = k)) +
  geom_line(size = 1.1) +
  geom_point(size = 2) +
  facet_grid(alpha ~ power, labeller = label_both, scales = "free_y") +
  # scale_color_grey(
  #   start = 0.1, end = 0.9,  # control how dark/light the grays are
  #   name = "k"
  # ) +
  scale_color_viridis_d(option = "turbo", name = "k") +
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
    legend.box = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 1))

###

# MINIMUM

pow_seq <- .8 #c(0.8, 0.9)
alpha_seq <- .05 #c(0.05, 0.01, 0.001)
k_seq <- 10 #c(3:10, 15, 20, 30)
decr_seq <- .2 #seq(0.1, 0.9, 0.1)
n_seq <- numeric(length(decr_seq))
n_master <- list()
for (pow in pow_seq) {
  for (alpha in alpha_seq) {
    cat("Power = ", pow, " - alpha = ", alpha, "\n", sep = "")
    n_all <- maxmin_multinom_size(k_seq, decr_seq, power = pow, alpha = alpha,
                                  type = "min", method = "Corrado", verbose = TRUE)
    n_master[[paste0("power = ", pow, " - alpha = ", alpha)]] <- n_all
  }
}

# debug(find_k_gamma)
# debug(maxmin_multinom_size)

save(n_master, file = "/Users/Sergio/dev/XOMultinom/data/sample_size_MIN.RData")

### TMP ###

# pow_seq <- 0.8
# alpha_seq <-0.05
# k_seq <- 40
# incr_seq <- 0.9
# probs_H0 <- rep(1/k_seq, k_seq)

# debugonce(find_k_gamma)
# nval <- 4
# find_k_gamma(probs_H0, nval, alpha_seq, type = "max", method = "Corrado")

# debugonce(maxmin_multinom_size)
# maxmin_multinom_size(k_seq, incr_seq, power = pow_seq, alpha = alpha_seq,
#                      type = "max", method = "Corrado", verbose = TRUE)
