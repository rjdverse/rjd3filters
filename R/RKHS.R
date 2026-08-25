#' Reproducing Kernel Hilbert Space (RKHS) Filters
#'
#' Estimation of a filter using Reproducing Kernel Hilbert Space (RKHS)
#' @inheritParams localpolynomials
#' @inheritParams mse
#' @param asymmetricCriterion the criteria used to compute the optimal
#'   bandwidth. If \code{"Undefined"}, \eqn{m+1} is used.
#' @param optimalbw boolean indicating if the bandwidth should be choosen by
#'   optimisation (between \code{optimal.minBandwidth} and
#'   \code{optimal.minBandwidth} using the criteria \code{asymmetricCriterion}).
#'   If \code{optimalbw = FALSE} then the bandwidth specified in
#'   \code{bandwidth} will be  used.
#' @param optimal.minBandwidth,optimal.maxBandwidth the range used for the
#'   optimal bandwidth selection.
#' @param bandwidth the bandwidth to use if \code{optimalbw = FALSE}.
#' @references Dagum, Estela Bee and Silvia Bianconcini (2008). “The Henderson
#'   Smoother in Reproducing Kernel Hilbert Space”. In: Journal of Business &
#'   Economic Statistics 26, pp. 536–545. URL:
#'   \url{https://ideas.repec.org/a/bes/jnlbes/v26y2008p536-545.html}.
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' rkhs <- rkhs_filter(horizon = 6, asymmetricCriterion = "Timeliness")
#' plot_coef(rkhs)
#' @returns A [finite_filters()] object.
#' @export
rkhs_filter <- function(
    horizon = 6,
    degree = 2,
    kernel = c(
        "BiWeight",
        "Henderson",
        "Epanechnikov",
        "Triangular",
        "Uniform",
        "TriWeight"
    ),
    asymmetricCriterion = c(
        "Timeliness",
        "FrequencyResponse",
        "Accuracy",
        "Smoothness",
        "Undefined"
    ),
    density = c("uniform", "rw"),
    passband = 2 * pi / 12,
    optimalbw = TRUE,
    optimal.minBandwidth = horizon,
    optimal.maxBandwidth = 3 * horizon,
    bandwidth = horizon + 1
) {
    kernel <- match.arg(
        tolower(kernel)[1],
        choices = c(
            "biweight",
            "henderson",
            "epanechnikov",
            "triangular",
            "uniform",
            "triweight"
        )
    )

    asymmetricCriterion <- switch(
        tolower(asymmetricCriterion[1]),
        timeliness = "Timeliness",
        frequencyresponse = "FrequencyResponse",
        accuracy = "Accuracy",
        smoothness = "Smoothness",
        undefined = "Undefined"
    )

    density <- match.arg(density)

    jrkhs_filter <-
        .jcall(
            "jdplus/filters/base/r/RKHSFilters",
            "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",
            "filters",
            as.integer(horizon),
            as.integer(degree),
            kernel,
            optimalbw,
            asymmetricCriterion,
            density == "rw",
            passband,
            bandwidth,
            optimal.minBandwidth,
            optimal.maxBandwidth
        )
    return(.jd2r_finitefilters(jrkhs_filter))
}
#' Optimization Function of Reproducing Kernel Hilbert Space (RKHS) Filters
#'
#' Export function used to compute the optimal bandwidth of Reproducing Kernel Hilbert Space (RKHS) filters
#' @inheritParams rkhs_filter
#' @inheritParams fst_filter
#' @returns A function that takes a bandwidth as input and returns the value of the optimization criterion.
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' graphics::plot(
#'     rkhs_optimization_fun(
#'         horizon = 6,
#'         leads = 0,
#'         degree = 3,
#'         asymmetricCriterion = "Timeliness"
#'     ),
#'     5.5,
#'     6 * 3,
#'     ylab = "Timeliness",
#'     main = "6X0 filter"
#' )
#' graphics::plot(
#'     rkhs_optimization_fun(
#'         horizon = 6,
#'         leads = 1,
#'         degree = 3,
#'         asymmetricCriterion = "Timeliness"
#'     ),
#'     5.5,
#'     6 * 3,
#'     ylab = "Timeliness",
#'     main = "6X1 filter"
#' )
#' graphics::plot(
#'     rkhs_optimization_fun(
#'         horizon = 6,
#'         leads = 2,
#'         degree = 3,
#'         asymmetricCriterion = "Timeliness"
#'     ),
#'     5.5,
#'     6 * 3,
#'     ylab = "Timeliness",
#'     main = "6X2 filter"
#' )
#' graphics::plot(
#'     rkhs_optimization_fun(
#'         horizon = 6,
#'         leads = 3,
#'         degree = 3,
#'         asymmetricCriterion = "Timeliness"
#'     ),
#'     5.5,
#'     6 * 3,
#'     ylab = "Timeliness",
#'     main = "6X3 filter"
#' )
#' graphics::plot(
#'     rkhs_optimization_fun(
#'         horizon = 6,
#'         leads = 4,
#'         degree = 3,
#'         asymmetricCriterion = "Timeliness"
#'     ),
#'     5.5,
#'     6 * 3,
#'     ylab = "Timeliness",
#'     main = "6X4 filter"
#' )
#' graphics::plot(
#'     rkhs_optimization_fun(
#'         horizon = 6,
#'         leads = 5,
#'         degree = 3,
#'         asymmetricCriterion = "Timeliness"
#'     ),
#'     5.5,
#'     6 * 3,
#'     ylab = "Timeliness",
#'     main = "6X5 filter"
#' )
#' @export
rkhs_optimization_fun <- function(
    horizon = 6,
    leads = 0,
    degree = 2,
    kernel = c(
        "Biweight",
        "Henderson",
        "Epanechnikov",
        "Triangular",
        "Uniform",
        "Triweight"
    ),
    asymmetricCriterion = c(
        "Timeliness",
        "FrequencyResponse",
        "Accuracy",
        "Smoothness"
    ),
    density = c("uniform", "rw"),
    passband = 2 * pi / 12
) {
    kernel <- match.arg(
        tolower(kernel)[1],
        choices = c(
            "biweight",
            "henderson",
            "epanechnikov",
            "triangular",
            "uniform",
            "triweight"
        )
    )
    asymmetricCriterion <- switch(
        tolower(asymmetricCriterion[1]),
        timeliness = "Timeliness",
        frequencyresponse = "FrequencyResponse",
        accuracy = "Accuracy",
        smoothness = "Smoothness",
        undefined = "Undefined"
    )
    density <- match.arg(density)
    jfun <-
        .jcall(
            "jdplus/filters/base/r/RKHSFilters",
            "Ljava/util/function/DoubleUnaryOperator;",
            "optimalCriteria",
            as.integer(horizon),
            as.integer(leads),
            as.integer(degree),
            kernel,
            asymmetricCriterion,
            density == "rw",
            passband
        )
    Vectorize(function(x) {
        .jcall(jfun, "D", "applyAsDouble", x)
    })
}
#' Optimal Bandwidth of Reproducing Kernel Hilbert Space (RKHS) Filters
#'
#' Function to compute the optimal bandwidths used in Reproducing Kernel Hilbert Space (RKHS) filters
#' @inheritParams rkhs_filter
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' rkhs_optimal_bw(asymmetricCriterion = "Timeliness")
#' rkhs_optimal_bw(asymmetricCriterion = "Timeliness", optimal.minBandwidth = 6.2)
#' @returns A vector of optimal bandwidths for each lead time.
#' @export
rkhs_optimal_bw <- function(
    horizon = 6,
    degree = 2,
    kernel = c(
        "Biweight",
        "Henderson",
        "Epanechnikov",
        "Triangular",
        "Uniform",
        "Triweight"
    ),
    asymmetricCriterion = c(
        "Timeliness",
        "FrequencyResponse",
        "Accuracy",
        "Smoothness"
    ),
    density = c("uniform", "rw"),
    passband = 2 * pi / 12,
    optimal.minBandwidth = horizon,
    optimal.maxBandwidth = 3 * horizon
) {
    kernel <- match.arg(
        tolower(kernel)[1],
        choices = c(
            "biweight",
            "henderson",
            "epanechnikov",
            "triangular",
            "uniform",
            "triweight"
        )
    )
    asymmetricCriterion <- switch(
        tolower(asymmetricCriterion[1]),
        timeliness = "Timeliness",
        frequencyresponse = "FrequencyResponse",
        accuracy = "Accuracy",
        smoothness = "Smoothness",
        undefined = "Undefined"
    )
    density <- match.arg(density)
    optimalBw <-
        .jcall(
            "jdplus/filters/base/r/RKHSFilters",
            "[D",
            "optimalBandwidth",
            as.integer(horizon),
            as.integer(degree),
            kernel,
            asymmetricCriterion,
            density == "rw",
            passband,
            optimal.minBandwidth,
            optimal.maxBandwidth
        )
    names(optimalBw) <- sprintf("q=%i", 0:(horizon - 1))
    optimalBw
}
#' Get RKHS kernel function
#' @inheritParams rkhs_filter
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' biweight <- rkhs_kernel(kernel = "Biweight")
#' triangular <- rkhs_kernel(kernel = "Triangular")
#' graphics::plot(biweight, -1, 1)
#' graphics::plot(triangular, -1, 1, add = TRUE, col = "orange")
#' @export A function that takes a numeric input and returns the value of the RKHS kernel.
#' @export
rkhs_kernel <- function(
    kernel = c(
        "Biweight",
        "Henderson",
        "Epanechnikov",
        "Triangular",
        "Uniform",
        "Triweight"
    ),
    degree = 2,
    horizon = 6
) {
    kernel <- match.arg(
        tolower(kernel)[1],
        choices = c(
            "biweight",
            "henderson",
            "epanechnikov",
            "triangular",
            "uniform",
            "triweight"
        )
    )
    kernel <- switch(
        tolower(kernel),
        biweight = "BiWeight",
        triweight = "TriWeight",
        uniform = "Uniform",
        triangular = "Triangular",
        epanechnikov = "Epanechnikov",
        henderson = "Henderson"
    )
    jfun <-
        .jcall(
            "jdplus/filters/base/r/RKHSFilters",
            "Ljava/util/function/DoubleUnaryOperator;",
            "kernel",
            kernel,
            as.integer(degree),
            as.integer(horizon)
        )
    Vectorize(function(x) {
        .jcall(jfun, "D", "applyAsDouble", x)
    })
}
