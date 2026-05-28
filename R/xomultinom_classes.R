# =============================================================================
# S3 classes and methods for XOMultinom
#
# Two classes are defined:
#
#   xomultinom_dist  -- returned by pmaxmultinom, dmaxmultinom, pminmultinom,
#                       dminmultinom, prangemultinom, drangemultinom,
#                       pJlargemultinom, dJlargemultinom.
#
#   xomultinom_size  -- returned by maxmin_multinom_size.
# =============================================================================


# =============================================================================
# CLASS: xomultinom_dist
# =============================================================================

#' Constructor for the \code{xomultinom_dist} class
#'
#' Internal constructor used by \code{\link{pmaxmultinom}},
#' \code{\link{dmaxmultinom}}, \code{\link{pminmultinom}},
#' \code{\link{dminmultinom}}, \code{\link{prangemultinom}}, and
#' \code{\link{drangemultinom}} to wrap their output in a structured S3 object.
#'
#' @param x Numeric vector of evaluation points.
#' @param values Numeric vector of probability (or log-probability) values of
#'   the same length as \code{x}.
#' @param stat Character string; one of \code{"max"}, \code{"min"},
#'   \code{"range"}, or \code{"J_largest"} indicating the order statistic.
#' @param type Character string; either \code{"pmf"} or \code{"cdf"}.
#' @param size Integer number of trials.
#' @param prob Numeric vector of (normalised) cell probabilities.
#' @param log Logical; \code{TRUE} if \code{values} are on the log scale.
#'
#' @return An object of class \code{xomultinom_dist}, which is a list with
#'   components \code{x}, \code{values}, \code{stat}, \code{type},
#'   \code{size}, \code{prob}, and \code{log}.
#'
#' @seealso \code{\link{print.xomultinom_dist}},
#'   \code{\link{summary.xomultinom_dist}},
#'   \code{\link{plot.xomultinom_dist}},
#'   \code{\link{autoplot.xomultinom_dist}}
#'
#' @noRd
#'
#' @keywords internal
new_xomultinom_dist <- function(x, values, stat, type, size, prob, log = FALSE) {
  stat <- match.arg(stat, c("max", "min", "range", "J_largest"))
  type <- match.arg(type, c("pmf", "cdf"))
  structure(
    list(
      x      = x,
      values = values,
      stat   = stat,
      type   = type,
      size   = size,
      prob   = prob,
      log    = log
    ),
    class = "xomultinom_dist"
  )
}


#' Print method for \code{xomultinom_dist} objects
#'
#' Displays a compact, human-readable table of evaluation points and the
#' corresponding exact probabilities (or log-probabilities) stored in an
#' \code{xomultinom_dist} object.
#'
#' @param x An object of class \code{xomultinom_dist}.
#' @param digits Integer number of significant digits for probabilities.
#'   Defaults to \code{4}.
#' @param max_rows Maximum number of rows to display when the support is large.
#'   If the number of evaluation points exceeds \code{max_rows}, the first and
#'   last \code{max_rows / 2} rows are shown with an ellipsis in between.
#'   Defaults to \code{20}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' k <- 5; n <- 40
#' obj <- dmaxmultinom(x = 0:n, size = n, prob = rep(1/k, k))
#' print(obj)
#'
#' @export
print.xomultinom_dist <- function(x, digits = 4, max_rows = 20, ...) {
  stat_label <- switch(x$stat,
    max       = "Maximum",
    min       = "Minimum",
    range     = "Range",
    J_largest = "J largest order stats"
  )
  type_label <- if (x$type == "pmf") "PMF" else "CDF"
  log_note   <- if (x$log) " [log scale]" else ""

  cat(sprintf(
    "Exact %s of the Multinomial %s%s\n",
    type_label, stat_label, log_note
  ))
  cat(sprintf(
    "  Trials (n) : %d    Categories (k) : %d\n",
    x$size, length(x$prob)
  ))

  if (length(unique(round(x$prob, 10))) == 1L) {
    cat(sprintf("  Probabilities  : equiprobable (p = %s)\n",
                format(x$prob[1], digits = digits)))
  } else {
    cat("  Probabilities  :", paste(round(x$prob, digits), collapse = ", "), "\n")
  }
  cat("\n")

  n_pts <- length(x$x)
  half  <- max_rows %/% 2L

  if (n_pts > max_rows) {
    idx_show <- c(seq_len(half), seq.int(n_pts - half + 1L, n_pts))
    x_show   <- x$x[idx_show]
    v_show   <- x$values[idx_show]
    cat(sprintf("  [displaying %d of %d evaluation points]\n\n", max_rows, n_pts))
  } else {
    x_show <- x$x
    v_show <- x$values
  }

  col_name <- if (x$type == "pmf") "P(X = x)" else "P(X <= x)"
  if (x$log) col_name <- paste0("log ", col_name)

  df <- data.frame(x = x_show, p = round(v_show, digits))
  names(df) <- c("x", col_name)

  if (n_pts > max_rows) {
    mid_row        <- data.frame(x = "...", p = "...")
    names(mid_row) <- names(df)
    df             <- rbind(df[seq_len(half), ],
                            mid_row,
                            df[seq.int(half + 1L, max_rows), ])
  }

  print(df, row.names = FALSE)
  invisible(x)
}


#' Summary method for \code{xomultinom_dist} objects
#'
#' Computes and displays descriptive statistics of the exact distribution
#' stored in an \code{xomultinom_dist} object, including the mean, median,
#' mode, standard deviation, effective support, and a central 95\% probability
#' interval.
#'
#' @param object An object of class \code{xomultinom_dist}.
#' @param digits Integer number of significant digits. Defaults to \code{4}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return Invisibly returns a named list with components \code{mean},
#'   \code{median}, \code{mode}, \code{sd}, \code{var}, \code{support},
#'   \code{q025}, and \code{q975}.
#'
#' @examples
#' k <- 5; n <- 40
#' obj <- pmaxmultinom(x = 0:n, size = n, prob = rep(1/k, k))
#' summary(obj)
#'
#' @export
summary.xomultinom_dist <- function(object, digits = 4, ...) {
  # Recover the PMF regardless of whether the object stores PMF or CDF
  if (object$type == "cdf") {
    vals <- if (object$log) exp(object$values) else object$values
    pmf  <- pmax(c(vals[1L], diff(vals)), 0)
  } else {
    pmf  <- if (object$log) exp(object$values) else object$values
  }
  pmf  <- pmf / sum(pmf)   # re-normalise to guard against floating-point drift
  xval <- object$x

  mean_val   <- sum(xval * pmf)
  var_val    <- sum((xval - mean_val)^2 * pmf)
  sd_val     <- sqrt(var_val)
  mode_val   <- xval[which.max(pmf)]

  cdf_vals   <- cumsum(pmf)
  median_val <- xval[which(cdf_vals >= 0.5)[1L]]
  q025       <- xval[which(cdf_vals >= 0.025)[1L]]
  q975       <- xval[which(cdf_vals >= 0.975)[1L]]
  supp       <- range(xval[pmf > .Machine$double.eps])

  stat_label <- switch(object$stat,
    max       = "Maximum",
    min       = "Minimum",
    range     = "Range",
    J_largest = "J largest order stats"
  )

  cat(sprintf("Exact distribution of the Multinomial %s\n", stat_label))
  cat(sprintf(
    "  n = %d    m = %d    %s\n\n",
    object$size, length(object$prob),
    if (length(unique(round(object$prob, 10))) == 1L) "equiprobable"
    else "non-equiprobable"
  ))
  cat(sprintf("  Mean              : %s\n", format(mean_val,   digits = digits)))
  cat(sprintf("  Median            : %s\n", format(median_val, digits = digits)))
  cat(sprintf("  Mode              : %s\n", format(mode_val,   digits = digits)))
  cat(sprintf("  Std deviation     : %s\n", format(sd_val,     digits = digits)))
  cat(sprintf("  Variance          : %s\n", format(var_val,    digits = digits)))
  cat(sprintf("  Effective support : [%s, %s]\n",
              format(supp[1L]), format(supp[2L])))
  cat(sprintf("  95%% interval      : [%s, %s]\n",
              format(q025), format(q975)))

  invisible(list(
    mean    = mean_val,
    median  = median_val,
    mode    = mode_val,
    sd      = sd_val,
    var     = var_val,
    support = supp,
    q025    = q025,
    q975    = q975
  ))
}


#' Plot method for \code{xomultinom_dist} objects
#'
#' Produces a base R plot of the exact distribution stored in an
#' \code{xomultinom_dist} object, compatible with \code{par(mfrow = ...)},
#' \code{layout()}, and all other base R multi-panel layout mechanisms.
#' PMFs are displayed as spike (needle) charts; CDFs are displayed as step
#' functions.  An optional normal approximation overlay can be added for
#' diagnostic comparison.
#'
#' @param x An object of class \code{xomultinom_dist}.
#' @param add_approx Logical; if \code{TRUE}, overlays the normal approximation
#'   to the distribution (mean and variance computed from the exact PMF).
#'   Defaults to \code{FALSE}.
#' @param col Character string; colour used for the exact distribution.
#'   Defaults to \code{"#2166ac"} (blue).
#' @param approx_col Character string; colour used for the approximation
#'   overlay when \code{add_approx = TRUE}. Defaults to \code{"#d6604d"}
#'   (red).
#' @param main Character string; plot title.  If \code{NULL} (default), a
#'   descriptive title is generated automatically.
#' @param xlab Character string; x-axis label. Defaults to \code{"x"}.
#' @param ylab Character string; y-axis label.  If \code{NULL} (default), an
#'   appropriate label is generated automatically.
#' @param ... Further graphical parameters passed to the underlying base R
#'   plotting functions.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @examples
#' k <- 5; n <- 40
#' obj_pmf <- dmaxmultinom(x = 0:n, size = n, prob = rep(1/k, k))
#' obj_cdf <- pmaxmultinom(x = 0:n, size = n, prob = rep(1/k, k))
#'
#' # Compatible with par(mfrow = ...)
#' op <- par(mfrow = c(1, 2))
#' plot(obj_pmf)
#' plot(obj_cdf)
#' par(op)
#'
#' @seealso \code{\link{autoplot.xomultinom_dist}} for a \code{ggplot2}-based
#'   alternative.
#'
#' @export
plot.xomultinom_dist <- function(x,
                                 add_approx = FALSE,
                                 col        = "#2166ac",
                                 approx_col = "#d6604d",
                                 main       = NULL,
                                 xlab       = "x",
                                 ylab       = NULL,
                                 ...) {
  stat_label <- switch(x$stat,
    max       = "Maximum",
    min       = "Minimum",
    range     = "Range",
    J_largest = "J largest order stats"
  )
  type_label <- if (x$type == "pmf") "PMF" else "CDF"

  if (is.null(main)) {
    main <- sprintf(
      "Multinomial %s - Exact %s  (n = %d, m = %d)",
      stat_label, type_label, x$size, length(x$prob)
    )
  }

  if (is.null(ylab)) {
    ylab <- if (x$type == "pmf") {
      if (x$log) "log P(X = x)" else "P(X = x)"
    } else {
      if (x$log) "log P(X <= x)" else "P(X <= x)"
    }
  }

  if (x$type == "pmf") {
    # Spike chart: plot empty frame first, then add segments and points
    plot(x$x, x$values,
         type = "n",
         main = main, xlab = xlab, ylab = ylab, ...)
    segments(x0  = x$x, y0 = 0,
             x1  = x$x, y1 = x$values,
             col = col, lwd = 1.5)
    points(x$x, x$values, pch = 19, cex = 0.8, col = col)
  } else {
    # Step function for the CDF
    plot(x$x, x$values,
         type = "s",
         col  = col, lwd = 1.5,
         main = main, xlab = xlab, ylab = ylab, ...)
  }

  # Optional normal approximation overlay
  if (add_approx && !x$log) {
    if (x$type == "cdf") {
      pmf_vals <- pmax(c(x$values[1L], diff(x$values)), 0)
    } else {
      pmf_vals <- x$values
    }
    pmf_vals <- pmf_vals / sum(pmf_vals)
    mu  <- sum(x$x * pmf_vals)
    sig <- sqrt(sum((x$x - mu)^2 * pmf_vals))

    x_seq    <- seq(min(x$x), max(x$x), length.out = 300)
    y_approx <- if (x$type == "pmf") {
      dnorm(x_seq, mean = mu, sd = sig)
    } else {
      pnorm(x_seq, mean = mu, sd = sig)
    }
    lines(x_seq, y_approx, col = approx_col, lwd = 1.5, lty = 2)
    legend("topleft",
           legend = c("Exact", "Normal approx."),
           col    = c(col, approx_col),
           lty    = c(1L, 2L),
           lwd    = 1.5,
           bty    = "n")
  }

  invisible(NULL)
}


#' ggplot2-based plot for \code{xomultinom_dist} objects
#'
#' Produces a \code{ggplot2} plot of the exact distribution stored in an
#' \code{xomultinom_dist} object.  PMFs are displayed as lollipop (spike)
#' charts; CDFs are displayed as step functions.  An optional normal
#' approximation overlay can be added for diagnostic comparison.
#'
#' For multi-panel layouts use \code{patchwork} or \code{gridExtra} to combine
#' multiple \code{autoplot()} outputs.  For base R \code{par(mfrow = ...)}
#' compatibility use \code{\link{plot.xomultinom_dist}} instead.
#'
#' @param object An object of class \code{xomultinom_dist}.
#' @param add_approx Logical; if \code{TRUE}, overlays the normal approximation
#'   to the distribution (mean and variance computed from the exact PMF).
#'   Defaults to \code{FALSE}.
#' @param colour Character string; colour used for the exact distribution.
#'   Defaults to \code{"#2166ac"} (blue).
#' @param approx_colour Character string; colour used for the approximation
#'   overlay when \code{add_approx = TRUE}. Defaults to \code{"#d6604d"}
#'   (red).
#' @param title Character string; plot title.  If \code{NULL} (default), a
#'   descriptive title is generated automatically.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @examples
#' k <- 5; n <- 40
#' obj <- dmaxmultinom(x = 0:n, size = n, prob = rep(1/k, k))
#' autoplot(obj)
#' autoplot(obj, add_approx = TRUE)
#'
#' @seealso \code{\link{plot.xomultinom_dist}} for a base R alternative
#'   compatible with \code{par(mfrow = ...)}.
#'
#' @importFrom ggplot2 autoplot
#' @export
autoplot.xomultinom_dist <- function(object,
                                     add_approx    = FALSE,
                                     colour        = "#2166ac",
                                     approx_colour = "#d6604d",
                                     title         = NULL,
                                     ...) {
  stat_label <- switch(object$stat,
    max       = "Maximum",
    min       = "Minimum",
    range     = "Range",
    J_largest = "J largest order stats"
  )
  type_label <- if (object$type == "pmf") "PMF" else "CDF"

  if (is.null(title)) {
    title <- sprintf(
      "Multinomial %s - Exact %s  (n = %d, m = %d)",
      stat_label, type_label, object$size, length(object$prob)
    )
  }

  y_label <- if (object$type == "pmf") {
    if (object$log) "log P(X = x)" else "P(X = x)"
  } else {
    if (object$log) "log P(X <= x)" else "P(X <= x)"
  }

  df <- data.frame(x = object$x, y = object$values)

  if (object$type == "pmf") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = .data$x, yend = 0),
        linewidth = 0.7, colour = colour
      ) +
      ggplot2::geom_point(size = 2, colour = colour)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_step(linewidth = 0.8, colour = colour)
  }

  # Optional normal approximation overlay
  if (add_approx && !object$log) {
    if (object$type == "cdf") {
      pmf_vals <- pmax(c(df$y[1L], diff(df$y)), 0)
    } else {
      pmf_vals <- df$y
    }
    pmf_vals <- pmf_vals / sum(pmf_vals)
    mu  <- sum(object$x * pmf_vals)
    sig <- sqrt(sum((object$x - mu)^2 * pmf_vals))

    x_seq    <- seq(min(object$x), max(object$x), length.out = 300)
    y_approx <- if (object$type == "pmf") {
      dnorm(x_seq, mean = mu, sd = sig)
    } else {
      pnorm(x_seq, mean = mu, sd = sig)
    }
    df_approx <- data.frame(x = x_seq, y = y_approx)
    p <- p +
      ggplot2::geom_line(
        data      = df_approx,
        ggplot2::aes(x = .data$x, y = .data$y),
        colour    = approx_colour,
        linewidth = 0.8,
        linetype  = "dashed"
      )
  }

  p <- p +
    ggplot2::labs(x = "x", y = y_label, title = title) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      axis.title = ggplot2::element_text(size = 11)
    )

  print(p)
  invisible(p)
}


#' Coerce an \code{xomultinom_dist} object to a data frame
#'
#' Converts the evaluation points and probability values stored in an
#' \code{xomultinom_dist} object into a tidy \code{data.frame} suitable for
#' further manipulation or export.
#'
#' @param x An object of class \code{xomultinom_dist}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return A \code{data.frame} with columns \code{x} (evaluation points) and
#'   either \code{pmf} or \code{cdf} (probability values). If the object was
#'   computed on the log scale the column is named \code{log_pmf} or
#'   \code{log_cdf} accordingly.
#'
#' @examples
#' k <- 5; n <- 40
#' obj <- dmaxmultinom(x = 0:n, size = n, prob = rep(1/k, k))
#' head(as.data.frame(obj))
#'
#' @export
as.data.frame.xomultinom_dist <- function(x, ...) {
  col_name      <- paste0(if (x$log) "log_" else "", x$type)
  df            <- data.frame(x = x$x, v = x$values)
  names(df)[2L] <- col_name
  df
}


# =============================================================================
# CLASS: xomultinom_size
# =============================================================================

#' Constructor for the \code{xomultinom_size} class
#'
#' Internal constructor used by \code{\link{maxmin_multinom_size}} to wrap its
#' output in a structured S3 object.
#'
#' @param sizes Named list of sample sizes as returned by
#'   \code{\link{maxmin_multinom_size}}.
#' @param m_seq Integer vector of numbers of categories used in the search.
#' @param change_seq Numeric vector of probability perturbations from
#'   equiprobability used in the search.
#' @param power Target power level.
#' @param alpha Significance level.
#' @param type Character string; either \code{"max"} or \code{"min"}.
#'
#' @return An object of class \code{xomultinom_size}.
#'
#' @seealso \code{\link{print.xomultinom_size}},
#'   \code{\link{summary.xomultinom_size}},
#'   \code{\link{plot.xomultinom_size}},
#'   \code{\link{autoplot.xomultinom_size}}
#'
#' @noRd
#'
#' @keywords internal
new_xomultinom_size <- function(sizes, m_seq, change_seq, power, alpha, type) {
  type <- match.arg(type, c("max", "min"))
  structure(
    list(
      sizes      = sizes,
      m_seq      = m_seq,
      change_seq = change_seq,
      power      = power,
      alpha      = alpha,
      type       = type
    ),
    class = "xomultinom_size"
  )
}


#' Print method for \code{xomultinom_size} objects
#'
#' Displays the required sample sizes as a formatted table, one block per
#' number of categories \eqn{m}.
#'
#' @param x An object of class \code{xomultinom_size}.
#' @param digits Integer number of decimal places for probability columns.
#'   Defaults to \code{4}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' sz <- maxmin_multinom_size(
#'   m_seq = c(5, 10), change_seq = c(0.10, 0.15, 0.20),
#'   power = 0.80, alpha = 0.05, type = "max"
#' )
#' print(sz)
#' }
#'
#' @export
print.xomultinom_size <- function(x, digits = 4, ...) {
  stat_label <- if (x$type == "max") "Maximum" else "Minimum"
  p_col      <- if (x$type == "max") "p_outlier" else "p_inlier"

  cat(sprintf(
    "Sample size requirements - Multinomial %s test\n",
    stat_label
  ))
  cat(sprintf(
    "  Target power : %.2f    Significance level : %.2f\n\n",
    x$power, x$alpha
  ))

  for (m in x$m_seq) {
    key         <- paste0("m = ", m)
    ns          <- x$sizes[[key]]
    cat(sprintf(" %s categories:\n", key))
    df          <- data.frame(as.numeric(names(ns)), as.integer(ns))
    names(df)   <- c(p_col, "n_required")
    df[[p_col]] <- round(df[[p_col]], digits)
    print(df, row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}


#' Summary method for \code{xomultinom_size} objects
#'
#' Prints a condensed overview of the required sample sizes across all
#' combinations of \eqn{m} and probability perturbations, reporting the
#' range of \eqn{n} for each \eqn{m}.
#'
#' @param object An object of class \code{xomultinom_size}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return Invisibly returns a named list where each element corresponds to a
#'   value of \code{m} and contains \code{n_min}, \code{n_max}, and
#'   \code{n_median}.
#'
#' @examples
#' \donttest{
#' sz <- maxmin_multinom_size(
#'   m_seq = c(5, 10), change_seq = c(0.10, 0.15, 0.20),
#'   power = 0.80, alpha = 0.05, type = "max"
#' )
#' summary(sz)
#' }
#'
#' @export
summary.xomultinom_size <- function(object, ...) {
  stat_label <- if (object$type == "max") "Maximum" else "Minimum"

  cat(sprintf(
    "Sample size summary - Multinomial %s test\n", stat_label
  ))
  cat(sprintf(
    "  Target power : %.2f    alpha : %.2f    m values : %s\n\n",
    object$power, object$alpha,
    paste(object$m_seq, collapse = ", ")
  ))

  out <- lapply(object$m_seq, function(m) {
    ns <- as.integer(object$sizes[[paste0("m = ", m)]])
    cat(sprintf(
      "  m = %2d :  n_min = %4d    n_max = %4d    n_median = %4d\n",
      m, min(ns), max(ns), as.integer(median(ns))
    ))
    list(n_min = min(ns), n_max = max(ns), n_median = as.integer(median(ns)))
  })
  names(out) <- paste0("m = ", object$m_seq)

  invisible(out)
}


#' Plot method for \code{xomultinom_size} objects
#'
#' Produces a base R line chart of the required sample size as a function of
#' the probability perturbation, with one line per value of \eqn{m} (number
#' of categories), compatible with \code{par(mfrow = ...)}, \code{layout()},
#' and all other base R multi-panel layout mechanisms.
#'
#' @param x An object of class \code{xomultinom_size}.
#' @param log_scale Logical; if \code{TRUE}, the \eqn{y}-axis (required
#'   \eqn{n}) is displayed on a \eqn{\log_{10}} scale. Useful when \eqn{n}
#'   varies over several orders of magnitude. Defaults to \code{FALSE}.
#' @param col Character vector of colours, one per value of \code{m_seq}.
#'   If \code{NULL} (default), colours are taken from the default R palette.
#' @param main Character string; plot title.  If \code{NULL} (default), a
#'   descriptive title is generated automatically.
#' @param xlab Character string; x-axis label. If \code{NULL} (default), an
#'   appropriate label is generated automatically.
#' @param ylab Character string; y-axis label. Defaults to
#'   \code{"Required n"}.
#' @param ... Further graphical parameters passed to the underlying base R
#'   plotting functions.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @examples
#' \donttest{
#' sz <- maxmin_multinom_size(
#'   m_seq = c(5, 10, 20), change_seq = seq(0.10, 0.30, by = 0.05),
#'   power = 0.80, alpha = 0.05, type = "max"
#' )
#'
#' # Compatible with par(mfrow = ...)
#' op <- par(mfrow = c(1, 2))
#' plot(sz)
#' plot(sz, log_scale = TRUE)
#' par(op)
#' }
#'
#' @seealso \code{\link{autoplot.xomultinom_size}} for a \code{ggplot2}-based
#'   alternative.
#'
#' @export
plot.xomultinom_size <- function(x,
                                 log_scale = FALSE,
                                 col       = NULL,
                                 main      = NULL,
                                 xlab      = NULL,
                                 ylab      = "Required n",
                                 ...) {
  stat_label <- if (x$type == "max") "Maximum" else "Minimum"

  if (is.null(main)) {
    main <- sprintf(
      "Multinomial %s test  (power = %.2f, alpha = %.2f)",
      stat_label, x$power, x$alpha
    )
  }
  if (is.null(xlab)) {
    xlab <- if (x$type == "max") "Outlier probability" else "Inlier probability"
  }
  if (is.null(col)) {
    col <- seq_along(x$m_seq)   # recycles through palette()
  }

  all_ns  <- unlist(lapply(x$m_seq, function(m)
    as.integer(x$sizes[[paste0("m = ", m)]])))
  log_arg <- if (log_scale) "y" else ""

  # Initialise empty plot with the correct axis ranges
  if (log_scale)
    ylab <- paste0(ylab, " (log-scale)")
  plot(range(x$change_seq), range(all_ns),
       type = "n",
       log  = log_arg,
       main = main, xlab = xlab, ylab = ylab, ...)

  for (i in seq_along(x$m_seq)) {
    m  <- x$m_seq[i]
    ns <- as.integer(x$sizes[[paste0("m = ", m)]])
    lines(x$change_seq, ns, col = col[i], lwd = 1.8, lty = i)
    points(x$change_seq, ns, col = col[i], pch = i, cex = 0.9)
  }

  legend("topright",
         legend = paste0("m = ", x$m_seq),
         col    = col,
         lty    = seq_along(x$m_seq),
         pch    = seq_along(x$m_seq),
         lwd    = 1.8,
         bty    = "n")

  invisible(NULL)
}


#' ggplot2-based plot for \code{xomultinom_size} objects
#'
#' Produces a \code{ggplot2} line chart of the required sample size as a
#' function of the probability perturbation, with one line per value of
#' \eqn{m} (number of categories).
#'
#' For multi-panel layouts use \code{patchwork} or \code{gridExtra} to combine
#' multiple \code{autoplot()} outputs.  For base R \code{par(mfrow = ...)}
#' compatibility use \code{\link{plot.xomultinom_size}} instead.
#'
#' @param object An object of class \code{xomultinom_size}.
#' @param log_scale Logical; if \code{TRUE}, the \eqn{y}-axis is displayed on
#'   a \eqn{\log_{10}} scale. Defaults to \code{FALSE}.
#' @param title Character string; plot title.  If \code{NULL} (default), a
#'   descriptive title is generated automatically.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' sz_max <- maxmin_multinom_size(
#'   m_seq = c(5, 10, 20), change_seq = seq(0.10, 0.30, by = 0.05),
#'   power = 0.80, alpha = 0.05, type = "max"
#' )
#' autoplot(sz_max)
#' autoplot(sz_max, log_scale = TRUE)
#' }
#'
#' @seealso \code{\link{plot.xomultinom_size}} for a base R alternative
#'   compatible with \code{par(mfrow = ...)}.
#'
#' @importFrom ggplot2 autoplot
#' @export
autoplot.xomultinom_size <- function(object,
                                     log_scale = FALSE,
                                     title     = NULL,
                                     ...) {
  stat_label <- if (object$type == "max") "Maximum" else "Minimum"
  x_label    <- if (object$type == "max") "Outlier probability" else "Inlier probability"

  if (is.null(title)) {
    title <- sprintf(
      "Required sample size - Multinomial %s test  (power = %.2f, alpha = %.2f)",
      stat_label, object$power, object$alpha
    )
  }

  rows <- lapply(object$m_seq, function(m) {
    ns <- object$sizes[[paste0("m = ", m)]]
    data.frame(
      m                = factor(paste0("m = ", m)),
      change           = as.numeric(names(ns)),
      n                = as.integer(ns),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)

  p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x      = .data$change,
        y      = .data$n,
        colour = .data$m,
        group  = .data$m
      )
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(
      x      = x_label,
      y      = ifelse(log_scale, "Required n (log-scale)", "Required n"),
      colour = "Categories",
      title  = title
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 12),
      axis.title      = ggplot2::element_text(size = 11),
      legend.position = "right"
    )

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }

  print(p)
  invisible(p)
}


#' Coerce an \code{xomultinom_size} object to a data frame
#'
#' Converts the sample size results stored in an \code{xomultinom_size} object
#' into a single tidy \code{data.frame} with columns for \eqn{m}, the
#' probability perturbation, and the required sample size.
#'
#' @param x An object of class \code{xomultinom_size}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return A \code{data.frame} with columns \code{m} (integer number of
#'   categories), \code{change} (probability perturbation), and
#'   \code{n_required} (required sample size).
#'
#' @examples
#' \donttest{
#' sz <- maxmin_multinom_size(
#'   m_seq = c(5, 10), change_seq = c(0.10, 0.15, 0.20),
#'   power = 0.80, alpha = 0.05, type = "max"
#' )
#' as.data.frame(sz)
#' }
#'
#' @export
as.data.frame.xomultinom_size <- function(x, ...) {
  rows <- lapply(x$m_seq, function(m) {
    ns <- x$sizes[[paste0("m = ", m)]]
    data.frame(
      m                = m,
      change           = as.numeric(names(ns)),
      n_required       = as.integer(ns),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
