dmaxmultinom <- function(x, size, prob, log = FALSE, verbose = FALSE, method = "Rcpp", parallel = "no",
  threads = 1, tol = 1e-10) {
  env <- environment()
	x_all <- numeric(2*length(x))
	x_all[seq(2, length(x_all), 2)] <- x
	x_all[seq(1, length(x_all), 2)] <- x - 1
	cdfmax_all <- pmaxmultinom(x_all, size, prob, log = FALSE, verbose = verbose, env = env, method = method,
    parallel = parallel, threads = threads, tol = tol)
  cdfmax_x <- cdfmax_all[seq(2, length(x_all), 2)]
  cdfmax_m1 <- cdfmax_all[seq(1, length(x_all), 2)]
  r <- cdfmax_x - cdfmax_m1
  
  if (!log) 
    r
  else log(r)
}
