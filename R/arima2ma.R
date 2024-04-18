#' Get Moving Averages from ARIMA model
#'
#' @param x the object.
#' @param ... unused parameters
#' @examples
#' fit <- stats::arima(log10(AirPassengers), c(0, 1, 1),
#' seasonal = list(order = c(0, 1, 1), period = 12))
#' get_moving_average(fit)
#' @export
get_moving_average <- function(x, ...) {
  UseMethod("get_moving_average", x)
}
#' @importFrom stats coefficients
#' @export
get_moving_average.Arima <- function(x, ...){
  arima_mod <- x$arma
  ar <- arima_mod[1]
  ma <- arima_mod[2]
  sar <- arima_mod[3]
  sma <- arima_mod[4]
  period <- arima_mod[5]
  diff <- arima_mod[6]
  sdiff <- arima_mod[7]

  ar_mm <- ma_mm <- sar_mm <-
    sma_mm <- moving_average(1, lags = 0)
  coef <- coefficients(x)
  if (ar > 0) {
    ar_mm <- moving_average(coef[sprintf("ar%i", seq(ar, 1))],
                            lags = - ar)
    ar_mm <- 1 - ar_mm
  }
  if (sar > 0) {
    sar_mm <- moving_average(coef[sprintf("sar%i", seq(sar, 1))],
                            lags = - sar)
    sar_mm <- to_seasonal(sar_mm, period)
    sar_mm <- 1 - sar_mm
  }
  if (ma > 0) {
    ma_mm <- moving_average(coef[sprintf("ma%i", seq(ma, 1))],
                            lags = - ma)
    ma_mm <- 1 + ma_mm
  }
  if (sma > 0) {
    sma_mm <- moving_average(coef[sprintf("sma%i", seq(sma, 1))],
                             lags = - sma)
    sma_mm <- to_seasonal(sma_mm, period)
    sma_mm <- 1 + sma_mm
  }

  # if (mean) {
  #   mean_mm <- x$regression.coefficients["Mean",1]
  # } else {
  #   mean_mm <- 0
  # }
  # mean_mm <- moving_average(mean_mm, 0)

  diff_mm <- (1 - moving_average(1, lags = -1)) ^ diff
  sdiff_mm <- (1 - moving_average(1, lags = -period)) ^ sdiff

  list(left = list(ar = ar_mm,
                 sar = sar_mm,
                 diff = diff_mm,
                 sdiff = sdiff_mm),
       right = list(ma = ma_mm,
                    sma = sma_mm))
}
#' @export
get_moving_average.regarima <- function(x, period = 12, ...){
  specif <- x$specification$arima$specification
  ar <- specif$arima.p
  ma <- specif$arima.q
  sar <- specif$arima.bp
  sma <- specif$arima.bq
  diff <- specif$arima.d
  sdiff <- specif$arima.bd
  mean <- x$model$spec_rslt$Mean

  ar_mm <- ma_mm <- sar_mm <-
    sma_mm <- moving_average(1, lags = 0)
  coef <- x$arima.coefficients[,1]

  if (ar > 0) {
    ar_mm <- moving_average(coef[sprintf("Phi(%i)", seq(ar, 1))],
                            lags = - ar)
    ar_mm <- 1 - ar_mm
  }
  if (sar > 0) {
    sar_mm <- moving_average(coef[sprintf("BPhi(%i)", seq(sar, 1))],
                             lags = - sar)
    sar_mm <- to_seasonal(sar_mm, period)
    sar_mm <- 1 - sar_mm
  }
  if (ma > 0) {
    ma_mm <- moving_average(coef[sprintf("Theta(%i)", seq(ma, 1))],
                            lags = - ma)
    ma_mm <- 1 - ma_mm
  }
  if (sma > 0) {
    sma_mm <- moving_average(coef[sprintf("BTheta(%i)", seq(sma, 1))],
                             lags = - sma)
    sma_mm <- to_seasonal(sma_mm, period)
    sma_mm <- 1 - sma_mm
  }
  if (mean) {
    mean_mm <- x$regression.coefficients["Mean",1]
  } else {
    mean_mm <- 0
  }
  mean_mm <- moving_average(mean_mm, 0)

  diff_mm <- (1 - moving_average(1, lags = -1)) ^ diff
  sdiff_mm <- (1 - moving_average(1, lags = -period)) ^ sdiff

  list(left = list(ar = ar_mm,
                   sar = sar_mm,
                   diff = diff_mm,
                   sdiff = sdiff_mm),
       right = list(ma = ma_mm,
                    sma = sma_mm))
}
#' @export
get_moving_average.SA <- function(x, period = 12, ...){
  get_moving_average(x$regarima, period = period, ...)
}
#' @export
get_moving_average.JD3_SARIMA_ESTIMATION <- function(x, period = 12, ...){
  ar <- x$phi
  ma <- x$theta
  sar <- x$bphi
  sma <- x$btheta
  diff <- x$d
  sdiff <- x$bd
  period <- x$period

  ar_mm <- ma_mm <- sar_mm <-
    sma_mm <- moving_average(1, lags = 0)

  if (! is.null(ar)) {
    ar_mm <- moving_average(rev(unlist(ar["value",])),
                            lags = - ncol(ar))
    ar_mm <- 1 - ar_mm
  }
  if (! is.null(sar)) {
    sar_mm <- moving_average(rev(unlist(sar["value",])),
                             lags = - ncol(sar))
    sar_mm <- to_seasonal(sar_mm, period)
    sar_mm <- 1 - sar_mm
  }
  if (! is.null(ma)) {
    ma_mm <- moving_average(rev(unlist(ma["value",])),
                            lags = - ncol(ma))
    ma_mm <- 1 - ma_mm
  }
  if (! is.null(sma)) {
    sma_mm <- moving_average(rev(unlist(sma["value",])),
                             lags = - ncol(sma))
    sma_mm <- to_seasonal(sma_mm, period)
    sma_mm <- 1 - sma_mm
  }

  diff_mm <- (1 - moving_average(1, lags = -1)) ^ diff
  sdiff_mm <- (1 - moving_average(1, lags = -period)) ^ sdiff

  list(left = list(ar = ar_mm,
                   sar = sar_mm,
                   diff = diff_mm,
                   sdiff = sdiff_mm),
       right = list(ma = ma_mm,
                    sma = sma_mm))
}
#' @export
get_moving_average.JD3_REGARIMA_OUTPUT <- function(x, ...){
  get_moving_average(x$result, ...)
}
#' @export
get_moving_average.JD3_REGARIMA_RSLTS <- function(x, ...){
  get_moving_average(x$description$arima, ...)
}
