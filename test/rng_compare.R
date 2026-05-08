# =============================================================================
# Benchmark: d-based vs p-based implementations for r functions
#
# Compares two implementations of rmaxmultinom (representative of all r
# functions) across a grid of (size, m, n_sim) values:
#
#   impl_d : current implementation -- uses dmaxmultinom() internally
#   impl_p : proposed implementation -- uses pmaxmultinom() + diff()
#
# Also included for reference:
#   impl_cdf: pure CDF inversion via runif() + .discrete_quantile()
#             (expected to be slowest for large n_sim)
#
# Requirements: XOMultinom, microbenchmark, ggplot2, dplyr
# =============================================================================

library(XOMultinom)
library(microbenchmark)
library(ggplot2)
library(dplyr)


# -----------------------------------------------------------------------------
# Helper: discrete quantile lookup (same as in qr_functions.R)
# -----------------------------------------------------------------------------
.discrete_quantile <- function(p, cdf, supp) {
  cdf[length(cdf)] <- 1
  vapply(p, function(pi) {
    idx <- which(cdf >= pi)
    if (length(idx) == 0L) supp[length(supp)] else supp[idx[1L]]
  }, integer(1L))
}


# -----------------------------------------------------------------------------
# Three implementations of rmaxmultinom
# -----------------------------------------------------------------------------

# Current implementation: uses dmaxmultinom (CDF -> diff inside C++ -> sample)
impl_d <- function(n_sim, size, prob) {
  supp <- 0L:size
  pmf  <- dmaxmultinom(x = supp, size = size, prob = prob,
                       log = FALSE, verbose = FALSE)$values
  pmf  <- pmax(pmf, 0)
  pmf  <- pmf / sum(pmf)
  sample(supp, size = n_sim, replace = TRUE, prob = pmf)
}

# Proposed implementation: uses pmaxmultinom + diff() in R -> sample
impl_p <- function(n_sim, size, prob) {
  supp <- 0L:size
  cdf  <- pmaxmultinom(x = supp, size = size, prob = prob,
                       log = FALSE, verbose = FALSE)$values
  pmf  <- c(cdf[1L], diff(cdf))
  pmf  <- pmax(pmf, 0)
  pmf  <- pmf / sum(pmf)
  sample(supp, size = n_sim, replace = TRUE, prob = pmf)
}

# Reference: pure CDF inversion (runif + quantile lookup; O(n_sim * size))
impl_cdf <- function(n_sim, size, prob) {
  supp <- 0L:size
  cdf  <- pmaxmultinom(x = supp, size = size, prob = prob,
                       log = FALSE, verbose = FALSE)$values
  .discrete_quantile(runif(n_sim), cdf, supp)
}


# -----------------------------------------------------------------------------
# Benchmark grid
# -----------------------------------------------------------------------------
grid <- expand.grid(
  size  = c(50L, 200L, 500L),   # number of trials
  m     = c(5L, 10L, 20L),      # number of categories
  n_sim = c(100L, 1000L, 5000L) # number of random draws
)

n_times <- 50L   # microbenchmark repetitions per cell

run_benchmark <- function(size, m, n_sim) {
  prob <- rep(1 / m, m)
  mb   <- microbenchmark(
    d   = impl_d(n_sim, size, prob),
    p   = impl_p(n_sim, size, prob),
    cdf = impl_cdf(n_sim, size, prob),
    times = n_times,
    unit  = "ms"
  )
  # Summarise: median time per implementation
  smry <- summary(mb)
  data.frame(
    size   = size,
    m      = m,
    n_sim  = n_sim,
    impl   = smry$expr,
    median = smry$median   # milliseconds
  )
}

cat("Running benchmarks ...\n")
results <- do.call(rbind, mapply(
  run_benchmark,
  size  = grid$size,
  m     = grid$m,
  n_sim = grid$n_sim,
  SIMPLIFY = FALSE
))
cat("Done.\n\n")


# -----------------------------------------------------------------------------
# Print summary table
# -----------------------------------------------------------------------------
tbl <- results |>
  tidyr::pivot_wider(names_from = impl, values_from = median) |>
  mutate(
    speedup_p   = round(d / p,   2),   # >1 means impl_p is faster
    speedup_cdf = round(d / cdf, 2)    # >1 means impl_cdf is faster
  ) |>
  arrange(size, m, n_sim)

cat("Median execution times (ms) and speedup relative to impl_d:\n\n")
print(tbl, digits = 3, row.names = FALSE)


# -----------------------------------------------------------------------------
# Plot 1: median time vs size, faceted by m, coloured by implementation
#         for a fixed n_sim
# -----------------------------------------------------------------------------
fixed_nsim <- 1000L

p1 <- results |>
  filter(n_sim == fixed_nsim) |>
  mutate(
    impl  = factor(impl,
                   levels = c("d", "p", "cdf"),
                   labels = c("impl_d  (dmaxmultinom)",
                              "impl_p  (pmaxmultinom + diff)",
                              "impl_cdf (pure inversion)")),
    m_lab = paste0("m = ", m)
  ) |>
  ggplot(aes(x = factor(size), y = median,
             colour = impl, group = impl)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ m_lab, nrow = 1) +
  scale_colour_manual(values = c("#2166ac", "#d6604d", "#4dac26")) +
  labs(
    x      = "size (number of trials)",
    y      = "Median time (ms)",
    colour = "Implementation",
    title  = sprintf("Random generation: median execution time vs size  (n_sim = %d)",
                     fixed_nsim),
    caption = "50 repetitions per cell; equiprobable multinomial"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92")
  )

print(p1)


# -----------------------------------------------------------------------------
# Plot 2: median time vs n_sim, faceted by m, for a fixed size
#         (shows where pure CDF inversion becomes expensive)
# -----------------------------------------------------------------------------
fixed_size <- 200L

p2 <- results |>
  filter(size == fixed_size) |>
  mutate(
    impl  = factor(impl,
                   levels = c("d", "p", "cdf"),
                   labels = c("impl_d  (dmaxmultinom)",
                              "impl_p  (pmaxmultinom + diff)",
                              "impl_cdf (pure inversion)")),
    m_lab = paste0("m = ", m)
  ) |>
  ggplot(aes(x = factor(n_sim), y = median,
             colour = impl, group = impl)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ m_lab, nrow = 1) +
  scale_colour_manual(values = c("#2166ac", "#d6604d", "#4dac26")) +
  labs(
    x      = "n_sim (number of random draws)",
    y      = "Median time (ms)",
    colour = "Implementation",
    title  = sprintf("Random generation: median execution time vs n_sim  (size = %d)",
                     fixed_size),
    caption = "50 repetitions per cell; equiprobable multinomial"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92")
  )

print(p2)


# -----------------------------------------------------------------------------
# Plot 3: speedup of impl_p over impl_d, as a heatmap over (size, m)
#         for fixed n_sim
# -----------------------------------------------------------------------------
p3 <- tbl |>
  filter(n_sim == fixed_nsim) |>
  mutate(across(c(size, m), factor)) |>
  ggplot(aes(x = size, y = m, fill = speedup_p)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f\u00d7", speedup_p)),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_gradient2(
    low      = "#d6604d",
    mid      = "grey80",
    high     = "#2166ac",
    midpoint = 1,
    name     = "Speedup\n(impl_p / impl_d)"
  ) +
  labs(
    x       = "size (number of trials)",
    y       = "m (number of categories)",
    title   = sprintf("Speedup of impl_p over impl_d  (n_sim = %d)", fixed_nsim),
    caption = "Values > 1 (blue) indicate impl_p is faster"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

print(p3)
