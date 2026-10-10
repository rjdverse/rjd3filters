#' Get Moving Averages from ARIMA model
#'
#' @param x the object.
#' @param ... unused parameters
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' fit <- stats::arima(log10(AirPassengers), c(0, 1, 1),
#' seasonal = list(order = c(0, 1, 1), period = 12))
#' get_moving_average(fit)
#'
#' @importFrom stats arima
#' @returns A [moving_average()] object.
#' @export
get_moving_average <- function(x, ...) {
    UseMethod("get_moving_average", x)
}

#' @noRd
#' @importFrom stats coefficients
#' @export
get_moving_average.Arima <- function(x, ...) {
    arima_mod <- x$arma
    order_ar <- arima_mod[1]
    order_ma <- arima_mod[2]
    order_sar <- arima_mod[3]
    order_sma <- arima_mod[4]
    period <- arima_mod[5]
    order_diff <- arima_mod[6]
    order_sdiff <- arima_mod[7]

    ar_mm <- ma_mm <- sar_mm <-
        sma_mm <- moving_average(1, lags = 0)
    mod_coef <- stats::coefficients(x)
    if (order_ar > 0) {
        ar_mm <- moving_average(
            mod_coef[sprintf("ar%i", seq(order_ar, 1))],
            lags = -order_ar
        )
        ar_mm <- 1 - ar_mm
    }
    if (order_sar > 0) {
        sar_mm <- moving_average(
            mod_coef[sprintf("sar%i", seq(order_sar, 1))],
            lags = -order_sar
        )
        sar_mm <- to_seasonal(sar_mm, period)
        sar_mm <- 1 - sar_mm
    }
    if (order_ma > 0) {
        ma_mm <- moving_average(
            mod_coef[sprintf("ma%i", seq(order_ma, 1))],
            lags = -order_ma
        )
        ma_mm <- 1 + ma_mm
    }
    if (order_sma > 0) {
        sma_mm <- moving_average(
            mod_coef[sprintf("sma%i", seq(order_sma, 1))],
            lags = -order_sma
        )
        sma_mm <- to_seasonal(sma_mm, period)
        sma_mm <- 1 + sma_mm
    }

    # if (mean) {
    #   mean_mm <- x$regression.coefficients["Mean",1]
    # } else {
    #   mean_mm <- 0
    # }
    # mean_mm <- moving_average(mean_mm, 0)

    diff_mm <- (1 - moving_average(1, lags = -1))^order_diff
    sdiff_mm <- (1 - moving_average(1, lags = -period))^order_sdiff

    list(
        left = list(ar = ar_mm, sar = sar_mm, diff = diff_mm, sdiff = sdiff_mm),
        right = list(ma = ma_mm, sma = sma_mm)
    )
}
#' @noRd
#' @export
get_moving_average.regarima <- function(x, period = 12, ...) {
    specif <- x$specification$arima$specification
    order_ar <- specif$arima.p
    order_ma <- specif$arima.q
    order_sar <- specif$arima.bp
    order_sma <- specif$arima.bq
    order_diff <- specif$arima.d
    order_sdiff <- specif$arima.bd
    order_mean <- x$model$spec_rslt$Mean

    ar_mm <- ma_mm <- sar_mm <-
        sma_mm <- moving_average(1, lags = 0)
    mod_coef <- x$arima.coefficients[, 1]

    if (order_ar > 0) {
        ar_mm <- moving_average(
            mod_coef[sprintf("Phi(%i)", seq(order_ar, 1))],
            lags = -order_ar
        )
        ar_mm <- 1 - ar_mm
    }
    if (order_sar > 0) {
        sar_mm <- moving_average(
            mod_coef[sprintf("BPhi(%i)", seq(order_sar, 1))],
            lags = -order_sar
        )
        sar_mm <- to_seasonal(sar_mm, period)
        sar_mm <- 1 - sar_mm
    }
    if (order_ma > 0) {
        ma_mm <- moving_average(
            mod_coef[sprintf("Theta(%i)", seq(order_ma, 1))],
            lags = -order_ma
        )
        ma_mm <- 1 - ma_mm
    }
    if (order_sma > 0) {
        sma_mm <- moving_average(
            mod_coef[sprintf("BTheta(%i)", seq(order_sma, 1))],
            lags = -order_sma
        )
        sma_mm <- to_seasonal(sma_mm, period)
        sma_mm <- 1 - sma_mm
    }
    if (order_mean) {
        mean_mm <- x$regression.coefficients["Mean", 1]
    } else {
        mean_mm <- 0
    }
    mean_mm <- moving_average(mean_mm, 0)

    diff_mm <- (1 - moving_average(1, lags = -1))^order_diff
    sdiff_mm <- (1 - moving_average(1, lags = -period))^order_sdiff

    list(
        left = list(ar = ar_mm, sar = sar_mm, diff = diff_mm, sdiff = sdiff_mm),
        right = list(ma = ma_mm, sma = sma_mm)
    )
}

#' @noRd
#' @export
get_moving_average.SA <- function(x, period = 12, ...) {
    get_moving_average(x$regarima, period = period, ...)
}
#' @noRd
#' @export
get_moving_average.JD3_SARIMA_ESTIMATION <- function(x, period = 12, ...) {
    order_ar <- x$phi
    order_ma <- x$theta
    order_sar <- x$bphi
    order_sma <- x$btheta
    order_diff <- x$d
    order_sdiff <- x$bd
    period <- x$period

    ar_mm <- ma_mm <- sar_mm <-
        sma_mm <- moving_average(1, lags = 0)

    if (!is.null(order_ar)) {
        ar_mm <- moving_average(
            rev(unlist(order_ar["value", ])),
            lags = -ncol(order_ar)
        )
        ar_mm <- 1 - ar_mm
    }
    if (!is.null(order_sar)) {
        sar_mm <- moving_average(
            rev(unlist(order_sar["value", ])),
            lags = -ncol(order_sar)
        )
        sar_mm <- to_seasonal(sar_mm, period)
        sar_mm <- 1 - sar_mm
    }
    if (!is.null(order_ma)) {
        ma_mm <- moving_average(
            rev(unlist(order_ma["value", ])),
            lags = -ncol(order_ma)
        )
        ma_mm <- 1 - ma_mm
    }
    if (!is.null(order_sma)) {
        sma_mm <- moving_average(
            rev(unlist(order_sma["value", ])),
            lags = -ncol(order_sma)
        )
        sma_mm <- to_seasonal(sma_mm, period)
        sma_mm <- 1 - sma_mm
    }

    diff_mm <- (1 - moving_average(1, lags = -1))^order_diff
    sdiff_mm <- (1 - moving_average(1, lags = -period))^order_sdiff

    list(
        left = list(ar = ar_mm, sar = sar_mm, diff = diff_mm, sdiff = sdiff_mm),
        right = list(ma = ma_mm, sma = sma_mm)
    )
}

#' @noRd
#' @export
get_moving_average.JD3_REGARIMA_OUTPUT <- function(x, ...) {
    get_moving_average(x$result, ...)
}

#' @noRd
#' @export
get_moving_average.JD3_REGARIMA_RSLTS <- function(x, ...) {
    get_moving_average(x$description$arima, ...)
}
