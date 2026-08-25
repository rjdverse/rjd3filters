#' Diagnostics and goodness of fit of filtered series
#'
#' Set of functions to compute diagnostics and goodness of fit of filtered series:
#' cross validation (`cv()`) and cross validate estimate (`cve()`),
#' leave-one-out cross validation estimate (`loocve`),
#' CP statistic (`cp()`) and Rice's T statistics (`rt()`).
#'
#' @param x input time series.
#' @param coef vector of coefficients or a moving-average ([moving_average()]).
#' @param ... other arguments passed to the function [moving_average()] to convert `coef` to a `"moving_average"` object.
#' @param var variance used to compute the CP statistic (`cp()`).
#'
#' @details Let \eqn{(\theta_i)_{-p\leq i \leq q}} be a moving average of length \eqn{p+q+1} used
#' to filter a time series \eqn{(y_i)_{1\leq i \leq n}}.
#' Let denote \eqn{\hat{\mu}_t} the filtered series computed at time \eqn{t} as:
#' \deqn{
#' \hat{\mu}_t = \sum_{i=-p}^q \theta_i y_{t+i}.
#' }
#'
#' The cross validation estimate (`cve()`) is defined as the time series \eqn{Y_t-\hat{\mu}_{-t}} where
#' \eqn{\hat{\mu}_{-t}} is the leave-one-out cross validation estimate (`loocve()`) defined as the filtered series
#' computed deleting the observation \eqn{t} and remaining all the other points.
#' The cross validation statistics (`cv()`) is defined as:
#' \deqn{
#' CV=\frac{1}{n-(p+q)}
#' \sum_{t=p+1}^{n-q} \left(y_t - \hat{\mu}_{-t}\right)^2.
#' }
#' In the case of filtering with a moving average, we can show that:
#' \deqn{
#' \hat{\mu}_{-t}= \frac{\hat{\mu}_t - \theta_0 y_t}{1-\theta_0}
#' }
#' and
#' \deqn{
#' CV=\frac{1}{n-(p+q)}
#' \sum_{t=p+1}^{n-q} \left(\frac{y_t - \hat{\mu}_{t}}{1-\theta_0}\right)^2.
#' }
#'
#' In the case of filtering with a moving average,
#' the CP estimate of risk (introduced by Mallows (1973); `cp()`) can be defined as:
#' \deqn{
#' CP=\frac{1}{\sigma^2}
#' \sum_{t=p+1}^{n-q} \left(y_t - \hat{\mu}_{t}\right)^2
#' -(n-(p+q))(1-2\theta_0).
#' }
#' The CP method requires an estimate of \eqn{\sigma^2} (`var` parameter).
#' The usual use of CP is to compare several different fits (for example different bandwidths):
#' one should use the same estimate of \eqn{\hat{\sigma}^2} for all fits (using for example [var_estimator()]).
#' The recommendation of Cleveland and Devlin (1988) is to compute \eqn{\hat{\sigma}^2}
#' from a fit at the smallest bandwidth under consideration,
#' at which one should be willing to assume that bias is negligible.
#'
#' The Rice's T statistic (`rt()`) is defined as:
#'\deqn{
#' \frac{1}{n-(p+q)}
#' \sum_{t=p+1}^{n-q}
#' \frac{
#'   \left(y_t - \hat{\mu}_{t}\right)^2
#' }{
#'   1-2\theta_0
#' }
#'}
#' @references
#' Loader, Clive. 1999.
#' Local regression and likelihood.
#' New York: Springer-Verlag.
#'
#' Mallows, C. L. (1973). Some comments on Cp.
#' Technometrics 15, 661– 675.
#'
#' Cleveland, W. S. and S. J. Devlin (1988).
#' Locally weighted regression: An approach to regression analysis by local fitting.
#' Journal of the American Statistical Association 83, 596–610.
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' x <- retailsa$DrinkingPlaces
#' h13 <- lp_filter(6)@sfilter
#' cv(x, h13)
#' cve(x, h13)
#' loocve(x, h13)
#' rt(x, h13)
#' cp(x, h13, var = var_estimator(x, h13))
#' @importFrom stats coef
#' @returns Numeric value for the statistic and `ts()` or `vector()` for the estimates.
#' @name diagnostics-fit
#' @rdname diagnostics-fit
#' @export
cve <- function(x, coef, ...) {
    coef <- moving_average(coef, ...)
    if (lower_bound(coef) > 0 || upper_bound(coef) < 0) {
        return(NA)
    }
    sc <- filter(x, coef)
    coef0 <- stats::coef(coef)["t"]
    (x - sc) / (1 - coef0)
}

#' @rdname diagnostics-fit
#' @export
cv <- function(x, coef, ...) {
    mean(cve(x, coef, ...)^2, na.rm = TRUE)
}

#' @importFrom stats coef
#' @rdname diagnostics-fit
#' @export
loocve <- function(x, coef, ...) {
    coef <- moving_average(coef, ...)
    if (lower_bound(coef) > 0 || upper_bound(coef) < 0) {
        return(NA)
    }
    sc <- filter(x, coef)
    coef0 <- stats::coef(coef)["t"]
    (sc - coef0 * x) / (1 - coef0)
}

#' @importFrom stats coef
#' @rdname diagnostics-fit
#' @export
rt <- function(x, coef, ...) {
    coef <- moving_average(coef, ...)
    if (lower_bound(coef) > 0 || upper_bound(coef) < 0) {
        return(NA)
    }
    sc <- filter(x, coef)
    coef0 <- stats::coef(coef)["t"]
    names(coef0) <- NULL
    mean((x - sc)^2, na.rm = TRUE) / (1 - 2 * coef0)
}

#' @importFrom stats coef
#' @rdname diagnostics-fit
#' @export
cp <- function(x, coef, var, ...) {
    mean(cve(x, coef, ...)^2, na.rm = TRUE)
    coef <- moving_average(coef, ...)
    if (lower_bound(coef) > 0 || upper_bound(coef) < 0) {
        return(NA)
    }
    sc <- filter(x, coef)
    coef0 <- stats::coef(coef)["t"]
    names(coef0) <- NULL
    mse <- (x - sc)^2
    nb_obs <- sum(!is.na(mse))
    (1 / var) * sum(mse, na.rm = TRUE) - nb_obs * (1 - 2 * coef0)
}

#' Variance Estimator
#'
#' @inheritParams diagnostics-fit
#'
#' @details
#' Let \eqn{(\theta_i)_{-p\leq i \leq q}} be a moving average of length \eqn{p+q+1} used
#' to filter a time series \eqn{(y_i)_{1\leq i \leq n}}.
#' It is equivalent to a local regression and the associated error variance \eqn{\sigma^2}
#' can be estimated using the normalized residual sum of squares, which can be simplified as:
#' \deqn{
#' \hat\sigma^2=\frac{1}{n-(p+q)}\sum_{t=p+1}^{n-q}
#' \frac{(y_t-\hat \mu_t)^2}{1-2w_0^2+\sum_{i=-p}^q w_i^2}
#' }
#' @references
#' Loader, Clive. 1999.
#' Local regression and likelihood.
#' New York: Springer-Verlag.
#' @seealso [df_var()].
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' x <- retailsa$DrinkingPlaces
#' h13 <- lp_filter(6)@sfilter
#' var_estimator(x, h13)
#'
#' @importFrom stats coefficients
#' @returns Numeric value with the statistic.
#' @export
var_estimator <- function(x, coef, ...) {
    coef <- moving_average(coef, ...)
    if (lower_bound(coef) > 0 || upper_bound(coef) < 0) {
        return(NA)
    }
    sc <- filter(x, coef)
    coef0 <- stats::coefficients(coef)["t"]
    sigma2 <- mean((x - sc)^2, na.rm = TRUE)
    sigma2 <- sigma2 / (1 - 2 * coef0 + sum(stats::coefficients(coef)^2))
    names(sigma2) <- NULL
    sigma2
}

#' Compute the degrees of freedom for the variance estimator
#'
#' @param n number of observations
#' @param coef moving average ([moving_average()]) used to filter the series.
#' @inheritParams confint_filter
#'
#' @seealso [var_estimator()].
#' @returns Numeric value with the statistic (the degree is not necessary an integer).
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' df_var(n = 100, coef = lp_filter(6)@sfilter)
#' @importFrom stats coefficients
#' @export
df_var <- function(n, coef, exact_df = FALSE) {
    value_coef <- stats::coefficients(coef)
    coef0 <- value_coef["t"]
    p <- abs(lower_bound(coef))
    f <- upper_bound(coef)
    df_num <- (n - (p + f)) * (1 - 2 * coef0 + sum(value_coef^2))
    names(df_num) <- NULL
    if (!exact_df) {
        return(df_num)
    } # Approximation of the degrees of freedom

    # Otherwise we compute the exact df more time consuming
    value_coef <- -value_coef
    value_coef["t"] <- 1 + value_coef["t"] # we already took the negative sign in the previous line
    mat_coefs <- do.call(
        cbind,
        lapply(0:(p + f), function(n_0) {
            c(
                rep(0, n_0),
                value_coef[seq(1, length.out = length(value_coef) - n_0)]
            )
        })
    )
    stats <- value_coef %*% mat_coefs
    stats <- stats^2
    stats[-1] <- stats[-1] * 2
    df_denum <- sum((n - (p + f) - seq(0, length.out = length(stats))) * stats)
    return(df_num^2 / df_denum)
}

#' Confidence intervals
#'
#' @inheritParams diagnostics-fit
#' @param coef moving-average ([moving_average()]) or finite filter
#'   ([finite_filters()]) used to filter the series.
#' @param coef_var moving-average ([moving_average()]) or finite filter
#'   ([finite_filters()]) used compute the variance (throw [var_estimator()]).
#' By default equal to `coef`.
#' @param level confidence level.
#' @param asymmetric_var if `asymmetric_var = TRUE` then the variance is
#'   estimated for each asymmetric filters instead of using the variance
#'   associated the symmetric estimates.
#' @param gaussian_distribution if `TRUE` use the normal distribution to
#'   compute the confidence interval, otherwise use the t-distribution.
#' @param exact_df if `TRUE` compute the exact degrees of freedom for the
#'   t-distribution (when `gaussian_distribution = FALSE`), otherwise uses an
#'   approximation.
#' @param ... other arguments passed to the function [moving_average()] to
#'   convert `coef` to a `"moving_average"` object.
#
#' @details
#' Let \eqn{(\theta_i)_{-p\leq i \leq q}} be a moving average of length
#' \eqn{p+q+1} used to filter a time series \eqn{(y_i)_{1\leq i \leq n}}.
#' Let denote \eqn{\hat{\mu}_t} the filtered series computed at time \eqn{t} as:
#' \deqn{
#' \hat{\mu}_t = \sum_{i=-p}^q \theta_i y_{t+i}.
#' }
#' If \eqn{\hat{\mu}_t} is unbiased, a approximate confidence for the true mean
#' is:
#' \deqn{
#' \left[\hat{\mu}_t - z_{1-\alpha/2} \hat{\sigma} \sqrt{\sum_{i=-p}^q\theta_i^2};
#' \hat{\mu}_t + z_{1-\alpha/2} \hat{\sigma} \sqrt{\sum_{i=-p}^q\theta_i^2}
#' \right],
#' }
#' where \eqn{z_{1-\alpha/2}} is the quantile \eqn{1-\alpha/2} of the standard
#' normal distribution.
#'
#' The estimate of the variance \eqn{\hat{\sigma}} is obtained using
#' [var_estimator()] with the parameter `coef_var`.
#' The assumption that \eqn{\hat{\mu}_t} is unbiased is rarely exactly true, so
#' variance estimates and confidence intervals are usually computed at small
#' bandwidths where bias is small.
#'
#' When `coef` (or `coef_var`) is a finite filter, the last points of the
#' confidence interval are computed using the corresponding asymmetric filters
#'
#' @references
#' Loader, Clive. 1999.
#' Local regression and likelihood.
#' New York: Springer-Verlag.
#'
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' x <- retailsa$DrinkingPlaces
#' coef <- lp_filter(6)
#' confint <- confint_filter(x, coef)
#' graphics::plot(confint, plot.type = "single",
#'      col = c("red", "black", "black"),
#'      lty = c(1, 2, 2))
#' @export
#' @returns A matrix or `mts` object with the filtered series and the lower and upper bounds of the confidence interval.
#' @importFrom stats frequency
#' @importFrom stats ts
#' @importFrom stats start
#' @importFrom stats qt
#' @importFrom stats qnorm
#' @importFrom stats end
#' @importFrom stats ts.union
#' @importFrom stats coefficients
confint_filter <- function(
    x,
    coef,
    coef_var = coef,
    level = 0.95,
    asymmetric_var = TRUE,
    gaussian_distribution = FALSE,
    exact_df = TRUE,
    ...
) {
    filtered <- filter(x, coef)
    confidence_bounds <- (1 - level) / 2
    confidence_bounds <- c(confidence_bounds, 1 - confidence_bounds)
    n <- length(filtered)
    if (is.moving_average(coef)) {
        corr_f <- sqrt(sum(stats::coefficients(coef)^2))
        if (gaussian_distribution) {
            confidence_quantiles <- matrix(
                stats::qnorm(confidence_bounds),
                ncol = 2
            )
        } else {
            confidence_quantiles <- matrix(
                stats::qt(
                    confidence_bounds,
                    df = df_var(n = n, coef = coef, exact_df = exact_df)
                ),
                ncol = 2
            )
        }
    } else if (is.finite_filters(coef)) {
        corr_f <- stats::ts(
            sqrt(sum(stats::coefficients(coef@sfilter)^2)),
            start = stats::start(filtered),
            end = stats::end(filtered),
            frequency = stats::frequency(filtered)
        )
        if (gaussian_distribution) {
            confidence_quantiles <- matrix(
                stats::qnorm(confidence_bounds),
                ncol = 2
            )
        } else {
            confidence_quantiles <- stats::ts(
                matrix(
                    stats::qt(
                        confidence_bounds,
                        df = df_var(
                            n = n,
                            coef = coef@sfilter,
                            exact_df = exact_df
                        )
                    ),
                    ncol = 2
                ),
                start = stats::start(filtered),
                end = stats::end(filtered),
                frequency = stats::frequency(filtered)
            )
        }
        lfilters <- coef@lfilters
        rfilters <- coef@rfilters
        for (i in seq_along(lfilters)) {
            corr_f[i] <- sqrt(sum(stats::coefficients(lfilters[[i]])^2))
            if (!gaussian_distribution) {
                confidence_quantiles[i, ] <- stats::qt(
                    confidence_bounds,
                    df = df_var(
                        n = n,
                        coef = lfilters[[i]],
                        exact_df = exact_df
                    )
                )
            }
        }
        for (i in seq_along(rfilters)) {
            corr_f[length(corr_f) - length(rfilters) + i] <-
                sqrt(sum(stats::coefficients(rfilters[[i]])^2))
            if (!gaussian_distribution) {
                confidence_quantiles[
                    length(time(confidence_quantiles)) - length(rfilters) + i,
                ] <-
                    stats::qt(
                        confidence_bounds,
                        df = df_var(
                            n = n,
                            coef = rfilters[[i]],
                            exact_df = exact_df
                        )
                    )
            }
        }
    }

    if (is.moving_average(coef_var)) {
        estimated_variance <- var_estimator(x, coef_var)
    } else if (is.finite_filters(coef_var)) {
        estimated_variance <- stats::ts(
            var_estimator(x, coef_var@sfilter),
            start = stats::start(filtered),
            end = stats::end(filtered),
            frequency = stats::frequency(filtered)
        )
        if (asymmetric_var) {
            lfilters <- coef_var@lfilters
            rfilters <- coef_var@rfilters
            for (i in seq_along(lfilters)) {
                estimated_variance[i] <- var_estimator(x, lfilters[[i]])
            }
            for (i in seq_along(rfilters)) {
                estimated_variance[
                    length(estimated_variance) - length(rfilters) + i
                ] <-
                    var_estimator(x, rfilters[[i]])
            }
        }
    }

    inf <- filtered +
        confidence_quantiles[, 1] * sqrt(estimated_variance) * corr_f
    sup <- filtered +
        confidence_quantiles[, 2] * sqrt(estimated_variance) * corr_f
    res <- base::cbind(filtered, inf, sup)
    colnames(res) <- c("filtered", sprintf("%.1f%%", confidence_bounds * 100))
    res
}
