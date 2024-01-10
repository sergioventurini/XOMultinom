round_exact <- function(x, digits = 0) {
  round(x + sign(x) * 1e-10, digits)
}
