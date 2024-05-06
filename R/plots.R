#' Plots filters properties
#'
#' Functions to plot the coefficients, the gain and the phase functions.
#'
#' @param x coefficients, gain or phase.
#' @param q q.
#' @param nxlab number of xlab.
#' @param ... other arguments to \code{matplot}.
#' @param add boolean indicating if the new plot is added to the previous one.
#' @param xlim vector containing x limits.
#' @param legend boolean indicating if the legend is printed.
#' @param legend.pos position of the legend.
#' @param normalized boolean indicatif if the phase function is normalized by the frequency.
#' @param zero_as_na boolean indicating if the trailing zero of the coefficients should be plotted (\code{FALSE}) or removed (\code{TRUE}).
#' @param n number of points used to plot the functions.
#'
#' @examples
#' filter <- lp_filter(6, endpoints = "DAF", kernel = "Henderson")
#' plot_coef(filter, q = c(0,3), legend = TRUE)
#' plot_gain(filter, q = c(0,3), legend = TRUE)
#' plot_phase(filter, q = c(0,3), legend = TRUE)
#' @name plot_filters
#' @rdname plot_filters
#' @importFrom MASS fractions
#' @export
plot_coef <- function(x, nxlab = 7, add = FALSE, ...){
  UseMethod("plot_coef", x)
}
#' @rdname plot_filters
#' @export
plot_coef.default <- function(x, nxlab = 7, add = FALSE,
                              zero_as_na = TRUE, q = 0, legend = FALSE,
                              legend.pos = "topright", ...){
  if (zero_as_na)
    x  <- apply(x,2, trailingzero_as_na)
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% colnames(x)]
  horizon <- (nrow(x)-1)/2
  if (length(col_to_plot) == 0) {
    if (!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "coefficient", xlim=c(-horizon, horizon), ylim=c(0, 1),
           ...)
      axis(1, at=seq(-horizon, horizon, by = 1), labels = rownames(x))
    }
    return(invisible(0))
  }
  matplot(seq(-horizon, horizon, by = 1),x[,col_to_plot],
          xaxt = "n", xlab = "", type = "o", pch = 20,
          ylab = "coefficient", add = add, ...)
  if (legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if (!add)
    axis(1, at=seq(-horizon, horizon, by = 1), labels = rownames(x))
}

#' @rdname plot_filters
#' @export
plot_coef.moving_average <- function(x, nxlab = 7, add = FALSE, ...){
  x_plot <- coef(x)
  matplot(seq(lower_bound(x), upper_bound(x), by = 1), x_plot,
          xaxt = "n", xlab = "", type = "o", pch = 20,
          ylab = "coefficient", add = add, ...)
  if (!add)
    axis(1, at=seq(lower_bound(x), upper_bound(x), by = 1), labels = names(x_plot))
}

#' @rdname plot_filters
#' @export
plot_coef.finite_filters <- function(x, nxlab = 7, add = FALSE,
                                     zero_as_na = TRUE, q = 0, legend = length(q) > 1,
                                     legend.pos = "topright", ...){
  plot_coef(x = as.matrix(x, zero_as_na = zero_as_na),
            nxlab = nxlab, add = add,
            zero_as_na = FALSE, q = q, legend = legend,
            legend.pos = legend.pos, ...)
}

#' @rdname plot_filters
#' @export
plot_gain <- function(x, nxlab = 7, add = FALSE,
                      xlim = c(0, pi), ...){
  UseMethod("plot_gain", x)
}
#' @rdname plot_filters
#' @export
plot_gain.moving_average<- function(x, nxlab = 7, add = FALSE,
                                    xlim = c(0, pi), ...){
  g <- get_properties_function(x, "Symmetric Gain")
  plot(g, type = "l",
       xaxt = "n", xlab = "",
       ylab = "gain", add = add, xlim = xlim, ...)
  if (!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}
#' @rdname plot_filters
#' @export
plot_gain.finite_filters <- function(x, nxlab = 7, add = FALSE,
                                xlim = c(0, pi), q = 0, legend = length(q) > 1,
                                legend.pos = "topright",
                                n = 101, ...){
  x_values <- seq.int(xlim[1], xlim[2], length.out = n)
  gsym <- get_properties_function(x, "Symmetric Gain")
  gasym <- get_properties_function(x, "Asymmetric Gain")
  all_g_f <- c(list(gsym), gasym)
  names(all_g_f)[1] <- sprintf("q=%i", upper_bound(x@sfilter))
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% names(all_g_f)]
  all_g_f <- all_g_f[col_to_plot]
  y_val <- sapply(all_g_f, function(f) f(x_values))
  if (length(col_to_plot) == 0) {
    if (!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "gain", xlim=xlim, ylim=c(0, 1),
           ...)
      x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
      axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
    }
    return(invisible(0))
  }
  matplot(x_values, y_val[, col_to_plot], type = "l",
          xaxt = "n", xlab = "",
          ylab = "gain", add = add, xlim = xlim, ...)

  if (legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if (!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}

#' @rdname plot_filters
#' @export
plot_phase <- function(x, nxlab = 7, add = FALSE,
                       xlim = c(0, pi), normalized = FALSE, ...){
  UseMethod("plot_phase", x)
}
#' @rdname plot_filters
#' @export
plot_phase.moving_average<- function(x, nxlab = 7, add = FALSE,
                                     xlim = c(0, pi), normalized = FALSE, ...){
  p <- get_properties_function(x, "Symmetric Phase")

  if (normalized) {
    p_plot <- function(x) {
      p(x)/x
      }
  } else {
    p_plot <- p
  }

  plot(p_plot, type = "l",
       xaxt = "n", xlab = "",
       ylab = "phase", add = add, xlim = xlim, ...)
  if (!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}


#' @rdname plot_filters
#' @export
plot_phase.finite_filters <- function(x, nxlab = 7, add = FALSE,
                                      xlim = c(0, pi), normalized = FALSE,
                                      q = 0, legend = length(q) > 1, legend.pos = "topright",
                                     n = 101, ...){
  x_values <- seq.int(xlim[1], xlim[2], length.out = n)
  psym <- get_properties_function(x, "Symmetric Phase")
  pasym <- get_properties_function(x, "Asymmetric Phase")
  all_p_f <- c(list(psym), pasym)
  names(all_p_f)[1] <- sprintf("q=%i", upper_bound(x@sfilter))
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% names(all_p_f)]
  all_p_f <- all_p_f[col_to_plot]
  y_val <- sapply(all_p_f, function(f) f(x_values))

  if (normalized){
    y_val[-1,] <- y_val[-1,] / x_values[-1]
  }
  if (length(col_to_plot) == 0) {
    if (!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "phase", xlim=xlim, ylim=c(0, 1),
           ...)
      x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
      axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
    }
    return(invisible(0))
  }
  matplot(x_values, y_val[, col_to_plot], type = "l",
          xaxt = "n", xlab = "",
          ylab = "phase", add = add, xlim = xlim, ...)

  if (legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if (!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}
xlabel <- function(x, symbol = "pi"){
  fracs <- strsplit(attr(MASS::fractions(x), "fracs"), "/")  # convert to fractions
  labels <- sapply(fracs, function(i)
    if (length(i) > 1) { paste(i[1], "*", symbol, "/", i[2]) }
    else { paste(i, "*", symbol) })
  labels <- sub("0 * pi", "0", labels, fixed = TRUE)
  labels <- sub("1 * pi", " pi", labels, fixed = TRUE)
  parse(text = labels)
}

trailingzero_as_na <- function(x){
  i <- length(x)
  while (x[i] == 0 && i > 0) {
    x[i] <- NA
    i <- i - 1
  }
  x
  # if (x[length(x)]==0)
  #   x [seq(from = tail(which(!sapply(x, function(y) isTRUE(all.equal(y,0)))),1)+1,
  #          to = length(x),
  #          by = 1)] <- NA
  # x
}
rm_leading_zero_or_na <- function(x){
  if (identical(x, 0))
    return(x)
  i <- 1
  remove_i <- NULL
  while ((is.na(x[i]) || (x[i] == 0)) && i <= length(x)) {
    remove_i <- c(i, remove_i)
    i <- i + 1
  }
  if (is.null(remove_i)){
    x
  } else{
    x[-remove_i]
  }
}
rm_trailing_zero_or_na <- function(x){
  if (identical(x, 0))
    return(x)
  i <- length(x)
  remove_i <- NULL
  while ((is.na(x[i]) || (x[i] == 0)) && i > 0) {
    remove_i <- c(i, remove_i)
    i <- i - 1
  }
  if (is.null(remove_i)){
    x
  } else{
    x[-remove_i]
  }
}
rm_trailing_zero <- function(x){
  if (identical(x, 0))
    return(x)
  i <- length(x)
  remove_i <- NULL
  while (isTRUE(all.equal(x[i], 0)) && i > 0) {
    remove_i <- c(i, remove_i)
    i <- i - 1
  }
  if (is.null(remove_i)){
    x
  } else{
    x[-remove_i]
  }
}
remove_bound_NA <- function(x) {
  if (all(is.na(x)))
    x
  i <- length(x)
  j <- 1
  remove_i_last <- remove_i_first <- NULL
  while (is.na(x[i]) && i > 0) {
    remove_i_last <- c(i, remove_i_last)
    i <- i - 1
  }
  while (is.na(x[j]) && i < length(x)) {
    remove_i_first <- c(j, remove_i_first)
    j <- j + 1
  }

  if (is.null(remove_i_first) & is.null(remove_i_last)){
    # list(data = x, leading = 0,
    #      trailing = 0)
  } else{
    x <- x[- c(remove_i_first, remove_i_last)]
  }

  list(data = x, leading = length(remove_i_first),
       trailing = length(remove_i_last))
}
