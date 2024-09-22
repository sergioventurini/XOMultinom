dminmultinom <- function(x, size, prob, log = FALSE, verbose = FALSE, method = "Rcpp") {
  if (all(diff(x) == 1)) {
  	if (method == "Rcpp") {
      cdfmin_all <- pminmultinom_C(x, size, prob, log = FALSE, verbose = verbose)
    } else if (method == "R") {
      cdfmin_all <- pminmultinom(x, size, prob, log = FALSE, verbose = verbose)
    } else {
      stop("method not available")
    }
    r <- c(cdfmin_all[1], diff(cdfmin_all))
  } else {
  	x_all <- numeric(2*length(x))
  	x_all[seq(2, length(x_all), 2)] <- x
  	x_all[seq(1, length(x_all), 2)] <- x - 1
    if (method == "Rcpp") {
    	cdfmin_all <- pminmultinom_C(x_all, size, prob, log = FALSE, verbose = verbose)
    } else if (method == "R") {
      cdfmin_all <- pminmultinom(x_all, size, prob, log = FALSE, verbose = verbose)
    } else {
      stop("method not available")
    }
    cdfmin_x <- cdfmin_all[seq(2, length(x_all), 2)]
    cdfmin_m1 <- cdfmin_all[seq(1, length(x_all), 2)]
    r <- cdfmin_x - cdfmin_m1
  }
  
  if (!log) 
    r
  else log(r)
}
