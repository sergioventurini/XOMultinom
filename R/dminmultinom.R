dminmultinom <- function(x, size, prob, log = FALSE, verbose = FALSE, method = "Rcpp", parallel = "no",
  threads = 1, tol = 1e-10) {
  env <- environment()
  if (all(diff(x) == 1)) {
    cdfmin_all <- pminmultinom(x, size, prob, log = FALSE, verbose = verbose, env = env, method = method,
      parallel = parallel, threads = threads, tol = tol)
    r <- c(cdfmin_all[1], diff(cdfmin_all))
  } else {
  	x_all <- numeric(2*length(x))
  	x_all[seq(2, length(x_all), 2)] <- x
  	x_all[seq(1, length(x_all), 2)] <- x - 1
  	cdfmin_all <- pminmultinom(x_all, size, prob, log = FALSE, verbose = verbose, env = env, method = method,
      parallel = parallel, threads = threads, tol = tol)
    cdfmin_x <- cdfmin_all[seq(2, length(x_all), 2)]
    cdfmin_m1 <- cdfmin_all[seq(1, length(x_all), 2)]
    r <- cdfmin_x - cdfmin_m1
  }
  
  if (!log) 
    r
  else log(r)
}
