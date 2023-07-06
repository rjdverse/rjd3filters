#' Compute IC-Ratio
#'
#' @param x input time series.
#' @param sc trend-cycle component.
#' @param mul boolean indicating if the decomposition is multiplicative or additive.
#'
#' @examples
#' x <- retailsa$AllOtherGenMerchandiseStores
#' sc <- henderson(x, length = 13, musgrave = FALSE)
#' ic_ratio(x, sc)
#'
#' @export
ic_ratio <- function(x, sc, mul = FALSE){
  remove_na <- is.na(x) | is.na(sc)
  x = as.numeric(x)[!remove_na]
  sc = as.numeric(sc)[!remove_na]
  result <- .jcall("jdplus/x11plus/base/r/X11Decomposition",
                   "D", "icratio",
                   x, sc, mul)
  result
}


# calcICR <- function(x, sc, length = 13){
#   si <- x - sc
#   gc <- calAbsMeanVariations(sc, 1)
#   gi <- calAbsMeanVariations(si, 1)
#   icr = gi/gc
#   freq = frequency(x)
#   if(freq == 4){
#     icr = icr * 3
#   } else if (freq == 2){
#     icr = icr * 6
#   }
#   return(icr)
# }
# calAbsMeanVariations <- function(x, nlags = 1, mul = FALSE){
#   mean = vector(mode = "numeric", length = nlags)
#   for (lag in 1:nlags){
#     x1 = x
#     x0 = lag(x, -lag)
#     d = x1 - x0
#     if(mul)
#       d = d/x0
#     mean[lag] = sum(abs(d),na.rm = TRUE)
#   }
#   mean
# }
#' X-11 Selection of Trend Filter
#'
#'
#' @param x I/C ratio [ic_ratio()] or a time series
#' @param ... further arguments passed to or from other methods.
#' @param freq frequency of the time series used to compute the I/C ratio.
#' @param length length of the Henderson filter used to compute the I/C ratio.
#'
#' @details The following procedure is used in X-11 to select the length of the trend filter:
#'
#' 1. Computes the I/C ratio, \eqn{icr} with an Henderson filter of length 13.
#'
#' 2. The length depends on the value or \eqn{icr}:
#'
#'    * if \eqn{icr < 1} then the selected length is 9 for monthly data and 5 otherwise;
#'    * if \eqn{1 \leq icr < 3.5} then the selected length is \eqn{freq + 1} where \eqn{freq} is the frequency of data (12 for monthly data, 4 for quarterly data...).
#'    * if \eqn{icr \geq 3.5} then the selected length is 23 for monthly data and 7 otherwise.
#' @examples
#' # example code
#' x <- retailsa$AllOtherGenMerchandiseStores
#' sc <- henderson(x, length = 13, musgrave = FALSE)
#' icr <- ic_ratio(x, sc)
#' select_trend_filter(icr, freq = 12)
#' # Because Henderson filter is used, this is equivalent to:
#' select_trend_filter(x)
#' @export
select_trend_filter <- function(x, ...){
  UseMethod("select_trend_filter")
}
#' @rdname select_trend_filter
#' @export
select_trend_filter.default <- function(x, ..., freq) {
  icr = x
  if (freq == 2) {
    return(c(icr = icr, length = 5));
  }
  if (icr >= 1 && icr < 3.5) {
    return(c(icr = icr, length = freq + 1));
  }
  if (icr < 1) {
    if (freq == 12) {
      c(icr = icr, length = 9)
    } else {
      c(icr = icr, length = 5)
    }
  } else {
    if(freq == 12){
      c(icr = icr, length = 23)
    }else{
      c(icr = icr, length = 7)
    }
  }
}
#' @rdname select_trend_filter
#' @export
select_trend_filter.ts <- function(x, ..., length = 13){
  select_trend_filter(ic_ratio(x, henderson(x, length = length, musgrave = FALSE)),
                        freq = frequency(x))
}
#' Cross-Validation
#'
#' Computes the cross-validation statistic of a moving-average on a time series
#'
#' @param x input time series.
#' @param coef vector of coefficients or a moving-average ([moving_average()]).
#' @param ... other arguments passed to the function [moving_average()] to convert `coef` to a `"moving_average"` object.
#'
#' @details The cross-validation statistics of a moving average \eqn{(\theta_i)_{-p\leq i \leq p}} on a time series \eqn{(y_i)_{1\leq i \leq n}} is:
#' \deqn{
#' \frac{1}{n-(p+q)}
#' \sum_{t=p+1}^{n-q}
#' \left(
#' \frac{y_t - \sum_{i=-p}^q \theta_i y_{t+i}}{
#' 1-\theta_0
#' }
#' \right)^2.
#' }
#' The function `cross_validation()` returns the time series \eqn{(x_t)}) with
#' \deqn{
#' x_t = \frac{y_t - \sum_{i=-p}^q \theta_i y_{t+i}}{
#' 1-\theta_0
#' }.
#' }
#'
#' @export
cross_validation <- function(x, coef, ...){
  coef <- moving_average(coef, ...)
  if (lower_bound(coef) > 0 || upper_bound(coef) < 0)
    return(NA)
  sc <- filter(x, coef)
  coef0 <- coef(coef)["t"]
  (x-sc)/(1-coef0)
}

#' Variance Estimator
#'
#' @inheritParams cross_validation
#'
#' @details
#' \deqn{
#' \hat\sigma^2=\frac{1}{n-2h}\sum_{t=h+1}^{n-h}\frac{(y_t-\hat \mu_t)^2}{1-2w_0^2+\sum w_i^2}
#' }
#'
#' @export
var_estimator <- function(x, coef, ...) {
  coef <- moving_average(coef, ...)
  if (lower_bound(coef) > 0 || upper_bound(coef) < 0)
    return(NA)
  sc <- filter(x, coef)
  coef0 <- coefficients(coef)["t"]
  sigma2 <-  mean((x - sc)^2,
                 na.rm = TRUE)
  sigma2 <- sigma2/(1- 2 * coef0 + sum(coef^2))
  names(sigma2) <- NULL
  sigma2
}
