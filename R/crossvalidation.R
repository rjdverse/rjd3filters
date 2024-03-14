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
  sigma2 <- sigma2/(1- 2 * coef0 + sum(coefficients(coef)^2))
  names(sigma2) <- NULL
  sigma2
}
