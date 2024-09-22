dmaxmultinom <- function(x, size, prob, log = FALSE, verbose = FALSE, method = "Rcpp") {
  if (all(diff(x) == 1)) {
    if (method == "Rcpp") {
    	cdfmax_all <- pmaxmultinom_C(x, size, prob, log = FALSE, verbose = verbose)
    } else if (method == "R") {
      cdfmax_all <- pmaxmultinom(x, size, prob, log = FALSE, verbose = verbose)
    } else {
      stop("method not available")
    }
    r <- c(cdfmax_all[1], diff(cdfmax_all))
  } else {
  	x_all <- numeric(2*length(x))
  	x_all[seq(2, length(x_all), 2)] <- x
  	x_all[seq(1, length(x_all), 2)] <- x - 1
    if (method == "Rcpp") {
    	cdfmax_all <- pmaxmultinom_C(x_all, size, prob, log = FALSE, verbose = verbose)
    } else if (method == "R") {
      cdfmax_all <- pmaxmultinom(x_all, size, prob, log = FALSE, verbose = verbose)
    } else {
      stop("method not available")
    }
    cdfmax_x <- cdfmax_all[seq(2, length(x_all), 2)]
    cdfmax_m1 <- cdfmax_all[seq(1, length(x_all), 2)]
    r <- cdfmax_x - cdfmax_m1
  }
  
  if (!log) 
    r
  else log(r)
}
