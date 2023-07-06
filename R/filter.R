#' Linear Filtering on a Time Series
#'
#' Applies linear filtering to a univariate time series or to each series separately of a multivariate time series using either a moving average (symmetric or asymmetric) or a combination of
#' symmetric moving average at the center and asymmetric moving averages at the bounds.
#'
#' @param x a univariate or multivariate time series.
#' @param coefs an object of class \code{"moving_average"} or the coefficients
#' of a moving average. In that case the argument \code{lags} must be provided.
#' @param coefs a \code{matrix} or a \code{list} that contains all the coefficients of the asymmetric and symmetric filters.
#'  (from the symmetric filter to the shortest). See details.
#'
#' @param remove_missing if `TRUE` (default) leading and trailing NA are removed before filtering.
#'
#' @details
#'
#'
#' The functions \code{filter} extends \code{\link[stats]{filter}} allowing to apply every kind of moving averages
#' (symmetric and asymmetric filters) or to apply aset multiple moving averages
#' to deal with the boundaries.
#'
#' Let \eqn{x_t} be the input time series to filter.
#'
#' * If `coef` is an object [moving_average()], of length \eqn{q}, the result \eqn{y} is equal at time \eqn{t} to:
#' \deqn{y[t] = x[t-lags] * coef[1] + x[t-lags+1] * coef[1] + ... + x[t-lags+q] * coef[q]}.
#' It extends the function \code{\link[stats]{filter}} that would add \code{NA} at the end of the time series.
#'
#' * If `coef` is a `matrix`, `list` or [finite_filters()] object,  at the center,
#' the symmetric moving average is used (first column/element of \code{coefs}).
#' At the boundaries, the last moving average of \code{coefs} is used to compute the filtered
#' time series \eqn{y[n]} (no future point known), the second to last to compute the filtered
#' time series \eqn{y[n-1]} (one future point known)...
#'
#' @examples
#' x <- retailsa$DrinkingPlaces
#'
#' lags <- 6
#' leads <- 2
#' fst_coef <- fst_filter(lags = lags, leads = leads, smoothness.weight = 0.3, timeliness.weight = 0.3)
#' lpp_coef <- lp_filter(horizon = lags, kernel = "Henderson", endpoints = "LC")
#'
#' fst_ma <- filter(x, fst_coef)
#' lpp_ma <- filter(x, lpp_coef[,"q=2"])
#'
#' plot(ts.union(x, fst_ma, lpp_ma), plot.type = "single", col = c("black","red","blue"))
#'
#' trend <- filter(x, lpp_coef)
#' # This is equivalent to:
#' trend <- localpolynomials(x, horizon = 6)
#' @export
filter <- function(x, coefs, remove_missing = TRUE){
  UseMethod("filter", x)
}
#' @export
filter.default <- function(x, coefs, remove_missing = TRUE){
  if (is.moving_average(coefs)) {
    filter_ma(x, coefs)
  } else {
    ff_ma(x, coefs = coefs, remove_missing = remove_missing)
  }
}
#' @export
filter.matrix <- function(x, coefs, remove_missing = TRUE){
  result <- x
  for (i in seq_len(ncol(x))){
    result[, i] <- filter(x[,i], coefs = coefs, remove_missing = remove_missing)
  }
  result
}

filter_ma <- function(x, coefs){
  # if (!is.moving_average(coefs)) {
  #   coefs <- moving_average(coefs, -abs(lags))
  # }
  lb = lower_bound(coefs)
  ub = upper_bound(coefs)

  if (length(x) <= length(coefs))
    return(x * NA)

  DataBlock = J("jdplus.toolkit.base.core.data.DataBlock")
  jx = DataBlock$of(as.numeric(x))
  out = DataBlock$of(as.numeric(rep(NA, length(x) - length(coefs)+1)))
  .ma2jd(coefs)$apply(jx,
                     out)
  result = out$toArray()
  result <- c(rep(NA, abs(min(lb, 0))),
              result,
              rep(NA, abs(max(ub, 0))))

  if (is.ts(x))
    result <- ts(result,start = start(x), frequency = frequency(x))
  result
}

ff_ma <- function(x, coefs, remove_missing = TRUE) {
  if (!inherits(coefs, "finite_filters")) {
    coefs <- finite_filters(coefs)
  }
  jffilters <- .finite_filters2jd(coefs)

  if (remove_missing) {
    data_clean <- remove_bound_NA(x)
    x2 <- data_clean$data
  } else {
    x2 <- x
  }

  jx <- .jcall("jdplus/toolkit/base/api/data/DoubleSeq",
               "Ljdplus/toolkit/base/api/data/DoubleSeq;",
               "of", as.numeric(x2))

  result <- .jcall("jdplus/toolkit/base/core/math/linearfilters/FilterUtility",
                   "Ljdplus/toolkit/base/api/data/DoubleSeq;", "filter",
                   jx,
                   jffilters$jsymf,
                   jffilters$jlasym,
                   jffilters$jrasym
  )

  result <- .jcall(result, "[D", "toArray")

  if (remove_missing){
    result = c(rep(NA, data_clean$leading), result,
               rep(NA, data_clean$trailing))
  }
  if(is.ts(x))
    result <- ts(result,start = start(x), frequency = frequency(x))
  result
}

.finite_filters2jd <- function(ff) {
  jsymf <- .ma2jd(ff@sfilter)
  rfilters <- ff@rfilters
  lfilters <- ff@lfilters
  if (length(rfilters) < upper_bound(ff@sfilter)) {
    # last points as NA
    rfilters <- c(rfilters,
                  rep(list(moving_average(NA, lags = 0)), upper_bound(ff@sfilter) - length(rfilters)))
  }
  if (length(lfilters) < -lower_bound(ff@sfilter)) {
    # first points as NA
    lfilters <- c(rep(list(moving_average(NA, lags = 0)), -lower_bound(ff@sfilter) - length(lfilters)),
                  lfilters)
  }
  jrasym <- lapply(rfilters, .ma2jd)
  jlasym <- rev(lapply(lfilters, .ma2jd))

  jsymf <- .jcast(jsymf,
                  "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter")
  if (length(jrasym) == 0) {
    jrasym <- .jnull("[Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;")
  } else {
    jrasym <- .jarray(jrasym,
                      "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter")
  }
  if (length(jlasym) == 0) {
    jlasym <- .jnull("[Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;")
  } else {
    jlasym <- .jarray(jlasym,
                      "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter")
  }
  list(jsymf = jsymf,
       jrasym = jrasym,
       jlasym = jlasym)
}
