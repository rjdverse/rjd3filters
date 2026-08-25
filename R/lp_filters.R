#' Apply Local Polynomials Filters
#'
#' @param x input time-series.
#' @param horizon horizon (bandwidth) of the symmetric filter.
#' @param degree degree of polynomial.
#' @param kernel kernel uses.
#' @param endpoints method for endpoints.
#' @param tweight timeliness weight.
#' @param passband passband threshold.
#' @param ic ic ratio.
#'
#' @returns A `vector` or `ts` object with the filtered time series.
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' x <- retailsa$AllOtherGenMerchandiseStores
#' trend <- localpolynomials(x, horizon = 6)
#' graphics::plot(x)
#' graphics::lines(trend, col = "red")
#' @references Proietti, Tommaso and Alessandra Luati (2008). “Real time estimation in local polynomial regression, with application to trend-cycle analysis”.
#' @seealso [lp_filter()].
#'
#' @importFrom stats frequency
#' @importFrom stats ts
#' @importFrom stats is.ts
#' @importFrom stats start
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @export
localpolynomials <- function(
    x,
    horizon = 6,
    degree = 3,
    kernel = c(
        "Henderson",
        "Uniform",
        "Biweight",
        "Trapezoidal",
        "Triweight",
        "Tricube",
        "Gaussian",
        "Triangular",
        "Parabolic"
    ),
    endpoints = c("LC", "QL", "CQ", "CC", "DAF"),
    ic = 3.5,
    tweight = 0,
    passband = pi / 12
) {
    if (2 * horizon < degree) {
        stop(
            "You need more observation (2 * horizon + 1) than variables (degree + 1) to estimate the filter."
        )
    }

    d <- 2 / (sqrt(pi) * ic)
    kernel <- match.arg(
        tolower(kernel),
        choices = c(
            "henderson",
            "uniform",
            "biweight",
            "trapezoidal",
            "triweight",
            "tricube",
            "gaussian",
            "triangular",
            "parabolic"
        )
    )
    kernel <- switch(
        kernel,
        henderson = "Henderson",
        uniform = "Uniform",
        biweight = "Biweight",
        trapezoidal = "Trapezoidal",
        triweight = "Triweight",
        tricube = "Tricube",
        gaussian = "Gaussian",
        triangular = "Triangular",
        parabolic = "Parabolic"
    )
    endpoints <- match.arg(endpoints)
    result <- .jcall(
        "jdplus/filters/base/r/LocalPolynomialFilters",
        "[D",
        "filter",
        as.numeric(x),
        as.integer(horizon),
        as.integer(degree),
        kernel,
        endpoints,
        d,
        tweight,
        passband
    )
    if (stats::is.ts(x)) {
        result <- stats::ts(
            result,
            start = stats::start(x),
            frequency = stats::frequency(x)
        )
    }
    result
}

#' Local Polynomials Filters
#'
#' @inheritParams localpolynomials
#' @details
#' * "LC": Linear-Constant filter
#' * "QL": Quadratic-Linear filter
#' * "CQ": Cubic-Quadratic filter
#' * "CC": Constant-Constant filter
#' * "DAF": Direct Asymmetric filter
#' * "CN": Cut and Normalized Filter
#'
#' @return a [finite_filters()] object.
#' @seealso [mmsre_filter()] [localpolynomials()].
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' henderson_f <- lp_filter(horizon = 6, kernel = "Henderson")
#' plot_coef(henderson_f)
#' @references Proietti, Tommaso and Alessandra Luati (2008). “Real time estimation in local polynomial regression, with application to trend-cycle analysis”.
#' @returns A [finite_filters()] object.
#' @export
lp_filter <- function(
    horizon = 6,
    degree = 3,
    kernel = c(
        "Henderson",
        "Uniform",
        "Biweight",
        "Trapezoidal",
        "Triweight",
        "Tricube",
        "Gaussian",
        "Triangular",
        "Parabolic"
    ),
    endpoints = c("LC", "QL", "CQ", "CC", "DAF", "CN"),
    ic = 3.5,
    tweight = 0,
    passband = pi / 12
) {
    if (2 * horizon < degree) {
        stop(
            "You need more observation (2 * horizon + 1) than variables (degree + 1) to estimate the filter."
        )
    }

    d <- 2 / (sqrt(pi) * ic)
    kernel <- match.arg(
        tolower(kernel),
        choices = c(
            "henderson",
            "uniform",
            "biweight",
            "trapezoidal",
            "triweight",
            "tricube",
            "gaussian",
            "triangular",
            "parabolic"
        )
    )
    kernel <- switch(
        kernel,
        henderson = "Henderson",
        uniform = "Uniform",
        biweight = "Biweight",
        trapezoidal = "Trapezoidal",
        triweight = "Triweight",
        tricube = "Tricube",
        gaussian = "Gaussian",
        triangular = "Triangular",
        parabolic = "Parabolic"
    )
    endpoints <- match.arg(endpoints)
    jprops <- .jcall(
        "jdplus/filters/base/r/LocalPolynomialFilters",
        "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",
        "filters",
        as.integer(horizon),
        as.integer(degree),
        kernel,
        endpoints,
        d,
        tweight,
        passband
    )

    return(.jd2r_finitefilters(jprops))
}
#' @keywords internal
coefficients_names <- function(lb, ub) {
    x <- sprintf("t%+i", seq(lb, ub))
    x <- sub("+0", "", x, fixed = TRUE)
    x
}
