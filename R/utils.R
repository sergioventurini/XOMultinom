round_exact <- function(x, digits = 0) {
  round(x + sign(x) * 1e-10, digits)
}

rdirichlet <- function (n, alpha) {
  n <- as.numeric(n)
  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow = 1)
  }
  if (prod(dim(alpha)) == 0) {
    stop2("alpha should be non-empty.")
  }
  if (isTRUE(any(alpha <= 0))) {
    stop2("alpha must be positive.")
  }
  if (n == 1) {
    n <- nrow(alpha)
  }
  if (n > nrow(alpha)) {
    alpha <- matrix(alpha, nrow = n, ncol = ncol(alpha), byrow = TRUE)
  }
  x <- matrix(rgamma(ncol(alpha) * n, alpha), ncol = ncol(alpha))
  x/rowSums(x)
}

list_2_array <- function(obj) {
  if (is.data.frame(obj)) {
    obj
  }
  obj_len <- length(names(obj))
  obj_nm_len <- nchar(names(obj)[1])
  dat <- data.frame()
  for (i in 1:obj_len) {
    dat <- rbind(dat, as.integer(unlist(strsplit(names(obj)[i], "_"))))
  }
  dat <- cbind(dat, unlist(obj, use.names = FALSE))
  rownames(dat) <- names(obj)
  colnames(dat) <- c(paste0("index_", 1:ceiling(obj_nm_len/2)), "prmax")
  x_seq_len <- max(dat[, 1]) + 1
  
  dat
}

row_match <- function(x, table, nomatch = NA) {
  if (inherits(table, "matrix"))
    table <- as.data.frame(table)
  if (is.null(dim(x)))
    x <- as.data.frame(matrix(x, nrow = 1))
  cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
  match(cx, ct, nomatch = nomatch)
}
