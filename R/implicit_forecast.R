#' Retrieve implicit forecasts corresponding to the asymmetric filters
#'
#' Function to retrieve the implicit forecasts corresponding to the asymmetric filters
#'
#' @details Let \eqn{h} be the bandwidth of the symmetric filter,
#' \eqn{v_{-h}, \ldots, v_h} the coefficients of the symmetric filter and
#' \eqn{w_{-h}^q, \ldots, w_h^q} the coefficients of the asymmetric filter used to estimate
#' the trend when \eqn{q} future values are known (with the convention \eqn{w_{q+1}^q=\ldots=w_h^q=0}).
#' Let denote \eqn{y_{-h},\ldots, y_0} the last \eqn{h} available values of the input times series.
#' The implicit forecast \eqn{y_{1}^*,\dots y_h^*} induced by \eqn{w^0,\dots w^{h-1}} are defined by:
#' \deqn{
#' \forall q\in\{0,...,h-1\}, \quad \sum_{i=-h}^0 v_iy_i + \sum_{i=1}^h v_iy_i^*
#' =\sum_{i=-h}^0 w_i^qy_i + \sum_{i=1}^h w_i^qy_i^*.
#' }
#' which is equivalent to
#' \deqn{
#' \forall q, \sum_{i=1}^h (v_i- w_i^q) y_i^*
#' =\sum_{i=-h}^0 (w_i^q-v_i)y_i.
#' }
#' Note that this is solved numerically: the solution isn't exact.
#' @inheritParams filter
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' x <- retailsa$AllOtherGenMerchandiseStores
#' ql <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "QL")
#' lc <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "LC")
#' f_ql <- implicit_forecasts(x, ql)
#' f_lc <- implicit_forecasts(x, lc)
#'
#' graphics::plot(window(x, start = 2007),
#'                xlim = c(2007,2012))
#' graphics::lines(
#'     x = stats::ts(
#'         c(utils::tail(x,1), f_ql),
#'         frequency = stats::frequency(x),
#'         start = stats::end(x)
#'     ),
#'     col = "red",
#'     lty = 2
#' )
#' graphics::lines(
#'     x = stats::ts(c(utils::tail(x,1), f_lc),
#'                   frequency = stats::frequency(x),
#'                   start = stats::end(x)
#'     ),
#'     col = "blue",
#'     lty = 2
#' )
#'
#' @importFrom stats frequency
#' @importFrom stats ts
#' @importFrom stats end
#' @importFrom utils tail
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @returns An object of the same class as `x` (`ts`, `mts`, `vector` or `matrix`) with the implicit forecasts.
#' @export
implicit_forecasts <- function(x, coefs) {
    UseMethod("implicit_forecasts", x)
}

#' @importFrom stats frequency
#' @importFrom stats ts
#' @importFrom stats deltat
#' @importFrom stats time
#' @importFrom stats is.ts
#' @importFrom utils tail
#' @noRd
#' @export
implicit_forecasts.default <- function(x, coefs) {
    if (!inherits(coefs, "finite_filters")) {
        coefs <- finite_filters(coefs)
    }
    jffilters <- .finite_filters2jd(coefs)

    jx <- .r2jd_doubleseq(utils::tail(x, abs(lower_bound(coefs@sfilter)) + 1))
    prev <- .jcall(
        "jdplus/toolkit/base/core/math/linearfilters/AsymmetricFiltersFactory",
        "[D",
        "implicitForecasts",
        jffilters$jsymf,
        jffilters$jrasym,
        jx
    )
    if (stats::is.ts(x)) {
        prev <- stats::ts(
            prev,
            frequency = stats::frequency(x),
            start = stats::time(x)[length(stats::time(x))] + stats::deltat(x)
        )
    }

    prev
}

#' @noRd
#' @export
implicit_forecasts.matrix <- function(x, coefs) {
    result <- do.call(
        cbind,
        lapply(
            seq_len(ncol(x)),
            function(i) implicit_forecasts(x[, i], coefs = coefs)
        )
    )
    colnames(result) <- colnames(x)
    result
}
