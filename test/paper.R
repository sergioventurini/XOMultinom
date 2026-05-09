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


### Section 5.1

library(ggplot2)
library(patchwork)

### Maximum test ###

pow_seq   <- c(0.8, 0.9)
alpha_seq <- c(0.05, 0.01, 0.001)
m_seq     <- c(3:10, 15, 20, 30)
incr_seq  <- seq(0.1, 0.9, 0.1)

params <- expand.grid(
  power = pow_seq,
  alpha = alpha_seq,
  KEEP.OUT.ATTRS = FALSE
)

n_master <- lapply(seq_len(nrow(params)), function(i) {
  pow <- params$power[i]  
  alpha <- params$alpha[i]

  cat("Power = ", pow, " - alpha = ", alpha, "\n", sep = "")
  
  res <- maxmin_multinom_size(
    m_seq, incr_seq,
    power = pow,
    alpha = alpha,
    n_max = 200,
    type = "max",
    verbose = TRUE,
    optmethod = "uniroot"
  )

  res$sizes
})

names(n_master) <- paste0(
  "power = ", params$power,
  " - alpha = ", params$alpha
)

df <- do.call(rbind, lapply(seq_along(n_master), function(i) {
  do.call(rbind, lapply(seq_along(m_seq), function(j) {
    data.frame(
      power = params$power[i],
      alpha = params$alpha[i],
      m     = m_seq[j],
      n     = n_master[[i]][[j]],
      incr  = incr_seq * 100
    )
  }))
}))

df$alpha <- factor(
  df$alpha,
  levels = sort(unique(df$alpha), decreasing = TRUE)
)

p_max <- ggplot(df, aes(x = incr, y = n, color = factor(m), group = m)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.5) +
  facet_grid(alpha ~ power, labeller = label_both, scales = "free_y") +
  scale_color_viridis_d(option = "turbo", name = "m") +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  labs(
    title = "Maximum test",
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
p_max <- p_max &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  
### Minimum test ###

pow_seq   <- c(0.8, 0.9)
alpha_seq <- c(0.05, 0.01, 0.001)
m_seq     <- 3:10
decr_seq  <- seq(0.2, 0.9, 0.1)

params <- expand.grid(
  power = pow_seq,
  alpha = alpha_seq,
  KEEP.OUT.ATTRS = FALSE
)

n_master <- lapply(seq_len(nrow(params)), function(i) {
  pow <- params$power[i]
  alpha <- params$alpha[i]

  cat("Power = ", pow, " - alpha = ", alpha, "\n", sep = "")

  res <- maxmin_multinom_size(
    m_seq, decr_seq,
    power = pow,
    alpha = alpha,
    n_max = 1000,
    type = "min",
    verbose = TRUE,
    optmethod = "uniroot"
  )

  res$sizes
})

names(n_master) <- paste0(
  "power = ", params$power,
  " - alpha = ", params$alpha
)

df <- do.call(rbind, lapply(seq_along(n_master), function(i) {
  do.call(rbind, lapply(seq_along(m_seq), function(j) {
    data.frame(
      power = params$power[i],
      alpha = params$alpha[i],
      m     = m_seq[j],
      n     = n_master[[i]][[j]],
      decr  = decr_seq * 100
    )
  }))
}))

df$alpha <- factor(
  df$alpha,
  levels = sort(unique(df$alpha), decreasing = TRUE)
)

p_min <- ggplot(df, aes(x = decr, y = n, color = factor(m), group = m)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.5) +
  facet_grid(alpha ~ power, labeller = label_both, scales = "free_y") +
  scale_color_viridis_d(option = "turbo", name = "m") +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  labs(
    title = "Minimum test",
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
p_min <- p_min &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

p_max + p_min


### Section 5.2

data(mainsail, package = "XOMultinom")

# Order data according to entry date
dat   <- mainsail[order(mainsail$entry_order), ]
score <- dat$halabi2014_lp_imputed
hist(
  score,
  col = "#4C78A8",
  border = "white",
  lwd = 1.2,
  freq = FALSE,
  xlab = "Score",
  ylab = "Density",
  main = "Halabi 2014 Risk Score for MAINSAIL data",
  xlim = c(-1, 1.6),
  ylim = c(0, 1),
  xaxt = "n",
  yaxt = "n"
)
axis(1, at = seq(-1, 16, 0.20), pos = 0)
axis(2, at = seq(0, 1, 0.10), pos = -1, las = 1)

# Overall settings
N     <- nrow(mainsail)  # MAINSAIL sample size
N0    <- 320L            # reference population size
n     <- 50L             # incoming batch size
m     <- 10L             # number of bins (deciles)
K     <- (N - N0)/n      # number of sequential tests
alpha <- 0.05            # nominal level
B     <- 1000L           # bootstrap replicates

set.seed(1406)

# Randomised-test parameters
p0      <- rep(1/m, m)
k_gamma <- find_k_gamma(probs = p0,
                        n = n,
                        alpha = alpha,
                        type = "max")
k_alpha <- k_gamma$k_alpha
gamma   <- k_gamma$gamma_prob

# Decile breaks (open outer boundaries)
make_breaks <- function(scores, m) {
  brks        <- quantile(scores, probs = seq(0, 1, by = 1/m), type = 1)
  brks[1]     <- -Inf
  brks[m + 1] <-  Inf
  unname(brks)
}

# Maximum bin count for a sample given breaks
max_count <- function(brks, samp_scores, m) {
  bins <- cut(samp_scores, breaks = brks, include.lowest = TRUE)
  max(tabulate(bins, nbins = m))
}

# Randomised test decision
rand_test <- function(obs_max, kappa, gamma) {
  if (obs_max >= kappa)      return(1L)
  if (obs_max == kappa - 1L) return(rbinom(1L, 1L, gamma))
  return(FALSE)
}

# Sequential procedure
results <- vector("list", K)
pop_idx <- seq_len(N0)
for (k in seq_len(K)) {
  samp_idx    <- (N0 + (k - 1L)*n + 1L):(N0 + k*n)
  ref_scores  <- score[pop_idx]
  samp_scores <- score[samp_idx]

  # observed test statistic and decision
  brks       <- make_breaks(ref_scores, m)
  counts_vec <- tabulate(cut(samp_scores, breaks = brks,
                             include.lowest = TRUE), nbins = m)
  obs_max    <- max(counts_vec)
  reject     <- rand_test(obs_max, k_alpha, gamma)

    # bootstrap sensitivity
  pi_hat <- mean(replicate(B, {
    b_brks <- make_breaks(sample(ref_scores, N0, replace = TRUE), m)
    rand_test(max_count(b_brks, samp_scores, m), k_alpha, gamma)
  }))

  results[[k]] <- list(
    pop   = paste0(pop_idx[1],  "\u2013", pop_idx[N0]),
    samp  = paste0(samp_idx[1], "\u2013", samp_idx[n]),
    counts = counts_vec,
    M     = obs_max,
    rej   = reject,
    pi    = pi_hat
  )

  # slide the window forward only on rejection
  if (reject)
    pop_idx <- c(pop_idx[-(1:n)], samp_idx)
}


panel_labels <- setNames(
  sapply(seq_along(results), function(k) {
    r <- results[[k]]
    sprintf("k = %d  \u00b7  cohort %s  \u00b7  M(%d) = %d  \u00b7  pi.hat = %.3f",
            k, r$samp, k, r$M, r$pi)
  }),
  paste0("k", seq_along(results))
)

df <- do.call(rbind, lapply(seq_along(results), function(k) {
  cnts <- results[[k]]$counts
  data.frame(
    test   = paste0("k", k),
    decile = factor(1:m),
    count  = cnts,
    is_max = cnts == max(cnts)
  )
}))
df$test  <- factor(df$test, levels = names(panel_labels))
df$label <- factor(panel_labels[as.character(df$test)], levels = panel_labels)

ref_lines <- data.frame(
  yintercept = c(k_alpha, n / m),
  linetype   = c("dashed", "dotted"),
  colour     = c("#C0392B", "#888780"),
  linewidth  = c(0.55, 0.45)
)

ref_labels <- c(
  sprintf("Critical value (k_alpha = %d)", k_alpha),
  sprintf("Expected under H0 (n/m = %d)", as.integer(n / m))
)

col_bar <- "#378ADD"
col_max <- "#EF9F27"

ggplot(df, aes(x = decile, y = count, fill = is_max)) +
  geom_col(width = 0.72, colour = NA) +
  geom_hline(
    data      = ref_lines,
    aes(yintercept = yintercept,
        linetype   = linetype,
        colour     = colour,
        linewidth  = linewidth),
    key_glyph = "path"
  ) +
  scale_fill_manual(
    values = c("FALSE" = col_bar, "TRUE" = col_max),
    labels = c("FALSE" = "Bin count", "TRUE" = "Maximum bin"),
    name   = NULL
  ) +
  scale_colour_identity(
    guide  = guide_legend(
      title        = NULL,
      override.aes = list(linetype  = c("dashed", "dotted"),
                          linewidth = c(0.55, 0.45))
    ),
    labels = ref_labels
  ) +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  scale_y_continuous(
    limits = c(0, k_alpha + 2),
    breaks = c(0, as.integer(n / m), k_alpha),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_x_discrete(breaks = as.character(seq_len(m))) +
  facet_wrap(~ label, ncol = 2,
             labeller = labeller(label = label_value)) +
  labs(x = "Decile bin", y = "Count") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
    strip.background   = element_rect(fill = "grey96", colour = "grey80",
                                      linewidth = 0.4),
    strip.text         = element_text(size = 8.5, margin = margin(4, 4, 4, 4)),
    legend.position    = "bottom",
    legend.key.width   = unit(1.6, "cm"),
    legend.text        = element_text(size = 8.5),
    legend.spacing.x   = unit(0.4, "cm"),
    axis.text          = element_text(size = 8),
    axis.title         = element_text(size = 9),
    plot.margin        = margin(6, 8, 4, 6)
  ) +
  guides(
    fill   = guide_legend(order = 1, override.aes = list(size = 4)),
    colour = guide_legend(order = 2)
  )
