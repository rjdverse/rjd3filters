#' Get properties of filters
#'
#' @param x a \code{"moving_average"} or \code{"finite_filters"} object.
#' @param component the component to extract.
#' @param ... unused other arguments.
#'
#' @examples
#' filter <- lp_filter(3, kernel = "Henderson")
#' sgain <- get_properties_function(filter, "Symmetric Gain")
#' plot(sgain, xlim= c(0, pi/12))
#' @export
get_properties_function <- function(x,
                                    component = c("Symmetric Gain",
                                                  "Symmetric Phase",
                                                  "Symmetric transfer",
                                                  "Asymmetric Gain",
                                                  "Asymmetric Phase",
                                                  "Asymmetric transfer"),
                                    ...) {
  UseMethod("get_properties_function", x)
}

get_gain_function <- function(x) {
  jgain <- .jcall(x, "Ljava/util/function/DoubleUnaryOperator;",
                  "gainFunction")
  Vectorize(function(x) {
    .jcall(jgain, "D", "applyAsDouble", x)
  })
}
get_phase_function <- function(x) {
  jphase <- .jcall(x, "Ljava/util/function/DoubleUnaryOperator;",
                   "phaseFunction")
  Vectorize(function(x) {
    .jcall(jphase, "D", "applyAsDouble", x)
  })
}
get_frequency_response_function <- function(x) {
  jfrf <- .jcall(x,
                 "Ljava/lang/Object;",
                 "frequencyResponseFunction")

  Vectorize(function(x) {
    res <- .jcall(jfrf, "Ljava/lang/Object;", "apply", x)

    complex(real = .jcall(res, "D", "getRe"),
            imaginary = .jcall(res, "D", "getIm"))
  })
}


#' @export
get_properties_function.moving_average <- function(x,
                                                   component = c("Symmetric Gain",
                                                                 "Symmetric Phase",
                                                                 "Symmetric transfer",
                                                                 "Asymmetric Gain",
                                                                 "Asymmetric Phase",
                                                                 "Asymmetric transfer"), ...) {
  x <- .ma2jd(x)
  component <- match.arg(component)
  switch(component,
         "Symmetric Gain" = {
           get_gain_function(x)
         },
         "Asymmetric Gain" = {
           get_gain_function(x)
         },
         "Symmetric Phase" = {
           get_phase_function(x)
         },
         "Asymmetric Phase" = {
           get_phase_function(x)
         },
         "Symmetric transfer" = {
           get_frequency_response_function(x)
         },
         "Asymmetric transfer" = {
           get_frequency_response_function(x)
         })
}
#' @export
get_properties_function.finite_filters <- function(x,
                                                   component = c("Symmetric Gain",
                                                                 "Symmetric Phase",
                                                                 "Symmetric transfer",
                                                                 "Asymmetric Gain",
                                                                 "Asymmetric Phase",
                                                                 "Asymmetric transfer"), ...) {
  component <- match.arg(component)
  if (any(grepl(pattern = "Symmetric", x = component, fixed = TRUE))) {
    get_properties_function(x@sfilter, component = component)
  } else {
    a_fun <- lapply(x@rfilters, get_properties_function, component = component)
    names(a_fun) <- sprintf("q=%i", seq(length(x@rfilters) - 1, 0))
    a_fun
  }
}

#' Compute quality criteria for asymmetric filters
#'
#' Function du compute a diagnostic matrix of quality criteria for asymmetric filters
#'
#' @param x Weights of the asymmetric filter (from -lags to m).
#' @param lags Lags of the filter (should be positive).
#' @param passband passband threshold.
#' @param sweights Weights of the symmetric filter (from 0 to lags or -lags to lags).
#' If missing, the criteria from the functions \code{\link{mse}} are not computed.
#' @param ... optional arguments to \code{\link{mse}}.
#'
#' @details For a moving average of coefficients \eqn{\theta=(\theta_i)_{-p\le i\le q}}
#' \code{diagnostic_matrix} returns a \code{list} with the following ten criteria:
#' \itemize{
#' \item{\code{b_c} } Constant bias (if \eqn{b_c=0}, \eqn{\theta} preserve constant trends)
#' \deqn{\sum_{i=-p}^q\theta_i - 1}
#' \item{\code{b_l} } Linear bias  (if \eqn{b_c=b_l=0}, \eqn{\theta} preserve constant trends)
#' \deqn{\sum_{i=-p}^q i \theta_i}
#' \item{\code{b_q} } Quadratic bias  (if \eqn{b_c=b_l=b_q=0}, \eqn{\theta} preserve quadratic trends)
#' \deqn{\sum_{i=-p}^q i^2 \theta_i}
#' \item{\code{F_g} } Fidelity criterium of Grun-Rehomme et al (2018)
#' \deqn{}
#' \item{\code{S_g} } Smoothness criterium of Grun-Rehomme et al (2018)
#' \item{\code{T_g} } Timeliness criterium of Grun-Rehomme et al (2018)
#' \item{\code{A_w} } Accuracy criterium of Wildi and McElroy (2019)
#' \item{\code{S_w} } Smoothness criterium of Wildi and McElroy (2019)
#' \item{\code{T_w} } Timeliness criterium of Wildi and McElroy (2019)
#' \item{\code{R_w} } Residual criterium of Wildi and McElroy (2019)
#' }
#'
#' @references Grun-Rehomme, Michel, Fabien Guggemos, and Dominique Ladiray (2018). “Asymmetric Moving Averages Minimizing Phase Shift”. In: Handbook on Seasonal Adjustment.
#'
#' Wildi, Marc and McElroy, Tucker (2019). “The trilemma between accuracy, timeliness and smoothness in real-time signal extraction”. In: International Journal of Forecasting 35.3, pp. 1072–1084.
#' @export
diagnostic_matrix <- function(x, lags, passband = pi/6,
                               sweights, ...) {
  if (!is.moving_average(x))
    x <- moving_average(x, lags = lags)

  results <- c(sum(x)-1, sum(coef(x) * seq(lower_bound(x), upper_bound(x), by = 1)),
               sum(coef(x) * seq(lower_bound(x), upper_bound(x), by = 1)^2),
               fst(x, lags, passband = passband))
  if (!missing(sweights)) {
    results <- c(results,
                 mse(x,
                     sweights,
                     passband = passband,
                     ...
                     )
    )
  } else {
    results <- c(results, rep(NA, 4))
  }
  names(results) <- c("b_c", "b_l", "b_q",
                      "F_g", "S_g", "T_g",
                      "A_w","S_w","T_w","R_w")
  results
}
