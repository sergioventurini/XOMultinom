drangemultinom <- function(x, size, prob, log = FALSE, verbose = FALSE, method = "Corrado", parallel = "no",
  threads = 1, tol = 1e-10) {
  env <- environment()
  if (all(diff(x) == 1)) {
  	cdfrange_all <- prangemultinom(x, size, prob, log = FALSE, verbose = verbose, env = env, method = method,
      parallel = parallel, threads = threads, tol = tol)
    r <- c(cdfrange_all[1], diff(cdfrange_all))
  } else {
  	x_all <- numeric(2*length(x))
  	x_all[seq(2, length(x_all), 2)] <- x
  	x_all[seq(1, length(x_all), 2)] <- x - 1
  	cdfrange_all <- prangemultinom(x_all, size, prob, log = FALSE, verbose = verbose, env = env, method = method,
      parallel = parallel, threads = threads, tol = tol)
    cdfrange_x <- cdfrange_all[seq(2, length(x_all), 2)]
    cdfrange_m1 <- cdfrange_all[seq(1, length(x_all), 2)]
    r <- cdfrange_x - cdfrange_m1
  }
  
  if (!log) 
    r
  else log(r)
}
