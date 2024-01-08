#' Reproducing Kernel Hilbert Space (RKHS) Filters
#'
#' Estimation of a filter using Reproducing Kernel Hilbert Space (RKHS)
#' @inheritParams localpolynomials
#' @inheritParams mse
#' @param asymmetricCriterion the criteria used to compute the optimal bandwidth. If \code{"Undefined"}, \eqn{m+1} is used.
#' @param optimalbw boolean indicating if the bandwith should be choosen by optimisation (between \code{optimal.minBandwidth} and
#' \code{optimal.minBandwidth} using the criteria \code{asymmetricCriterion}).
#' If \code{optimalbw = FALSE} then the bandwith specified in \code{bandwidth} will be  used.
#' @param optimal.minBandwidth,optimal.maxBandwidth the range used for the optimal bandwith selection.
#' @param bandwidth the bandwidth to use if \code{optimalbw = FALSE}.
#' @references Dagum, Estela Bee and Silvia Bianconcini (2008). “The Henderson Smoother in Reproducing Kernel Hilbert Space”. In: Journal of Business & Economic Statistics 26, pp. 536–545. URL: \url{https://ideas.repec.org/a/bes/jnlbes/v26y2008p536-545.html}.
#' @examples
#' rkhs <- rkhs_filter(horizon = 6, asymmetricCriterion = "Timeliness")
#' plot_coef(rkhs)
#' @return An object of class \code{"rkhs_filter"}, which is a list of 4 elements:\itemize{
#' \item{\code{"internal"}}{Java object used for internal functions}
#' \item{\code{"filters.coef"}}{The coefficients of the selected filter}
#' \item{\code{"filters.gain"}}{The gain function between 0 and pi (601 observations)}
#' \item{\code{"filters.phase"}}{The phase function between 0 and pi (601 observations)}
#' }
#' @export
rkhs_filter <- function(horizon = 6, degree = 2,
                        kernel = c("BiWeight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "TriWeight"),
                        asymmetricCriterion = c("Timeliness", "FrequencyResponse", "Accuracy", "Smoothness", "Undefined"),
                        density = c("uniform", "rw"),
                        passband = 2*pi/12,
                        optimalbw = TRUE,
                        optimal.minBandwidth = horizon,
                        optimal.maxBandwidth = 3*horizon,
                        bandwidth = horizon + 1){
  kernel <- match.arg(tolower(kernel)[1],
                   choices = c("biweight", "henderson", "epanechnikov", "triangular", "uniform",
                               "triweight"))
  # In next release of Java files, remove next lines
  kernel = switch(tolower(kernel),
                  "biweight" = "BiWeight",
                  "triweight" ="TriWeight",
                  "uniform" = "Uniform",
                  "triangular" = "Triangular",
                  "epanechnikov" = "Epanechnikov",
                  "henderson" = "Henderson"
  )

  asymmetricCriterion = switch(tolower(asymmetricCriterion[1]),
                               timeliness = "Timeliness",
                               frequencyresponse = "FrequencyResponse",
                               accuracy = "Accuracy",
                               smoothness = "Smoothness",
                               undefined = "Undefined")

  density = match.arg(density)

  jrkhs_filter =
    .jcall("jdplus/filters/base/r/RKHSFilters",
           "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",
           "filters",
    as.integer(horizon), as.integer(degree), kernel,
    optimalbw, asymmetricCriterion, density=="rw", passband,
    bandwidth, optimal.minBandwidth, optimal.maxBandwidth
  )
  return(.jd2r_finitefilters(jrkhs_filter))
}
#' Optimization Function of Reproducing Kernel Hilbert Space (RKHS) Filters
#'
#' Export function used to compute the optimal bandwidth of Reproducing Kernel Hilbert Space (RKHS) filters
#' @inheritParams rkhs_filter
#' @inheritParams fst_filter
#' @examples
#' plot(rkhs_optimization_fun(horizon = 6, leads = 0,degree = 3, asymmetricCriterion = "Timeliness"),
#'      5.5, 6*3, ylab = "Timeliness",
#'      main = "6X0 filter")
#' plot(rkhs_optimization_fun(horizon = 6, leads = 1,degree = 3, asymmetricCriterion = "Timeliness"),
#'      5.5, 6*3, ylab = "Timeliness",
#'      main = "6X1 filter")
#' plot(rkhs_optimization_fun(horizon = 6, leads = 2,degree = 3, asymmetricCriterion = "Timeliness"),
#'      5.5, 6*3, ylab = "Timeliness",
#'      main = "6X2 filter")
#' plot(rkhs_optimization_fun(horizon = 6, leads = 3,degree = 3, asymmetricCriterion = "Timeliness"),
#'      5.5, 6*3, ylab = "Timeliness",
#'      main = "6X3 filter")
#' plot(rkhs_optimization_fun(horizon = 6, leads = 4,degree = 3, asymmetricCriterion = "Timeliness"),
#'      5.5, 6*3, ylab = "Timeliness",
#'      main = "6X4 filter")
#' plot(rkhs_optimization_fun(horizon = 6, leads = 5,degree = 3, asymmetricCriterion = "Timeliness"),
#'      5.5, 6*3, ylab = "Timeliness",
#'      main = "6X5 filter")
#' @export
rkhs_optimization_fun <- function(horizon = 6, leads = 0,  degree = 2,
                        kernel = c("Biweight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "Triweight"),
                        asymmetricCriterion = c("Timeliness", "FrequencyResponse", "Accuracy", "Smoothness"),
                        density = c("uniform", "rw"),
                        passband = 2*pi/12){
  kernel <- match.arg(tolower(kernel)[1],
                   choices = c("biweight", "henderson", "epanechnikov", "triangular", "uniform",
                               "triweight"))
  asymmetricCriterion = switch(tolower(asymmetricCriterion[1]),
                               timeliness = "Timeliness",
                               frequencyresponse = "FrequencyResponse",
                               accuracy = "Accuracy",
                               smoothness = "Smoothness",
                               undefined = "Undefined")
  density = match.arg(density)
  optimalFunCriteria = J("jdplus/filters/base/r/RKHSFilters")$optimalCriteria(
    as.integer(horizon), as.integer(leads), as.integer(degree), kernel,
    asymmetricCriterion, density=="rw", passband
  )$applyAsDouble

  Vectorize(function(x){
    optimalFunCriteria(x)
  })
}
#' Optimal Bandwith of Reproducing Kernel Hilbert Space (RKHS) Filters
#'
#' Function to export the optimal bandwidths used in Reproducing Kernel Hilbert Space (RKHS) filters
#' @inheritParams rkhs_filter
#' @examples
#' rkhs_optimal_bw(asymmetricCriterion = "Timeliness")
#' rkhs_optimal_bw(asymmetricCriterion = "Timeliness", optimal.minBandwidth = 6.2)
#' @export
rkhs_optimal_bw <- function(horizon = 6,  degree = 2,
                           kernel = c("Biweight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "Triweight"),
                           asymmetricCriterion = c("Timeliness", "FrequencyResponse", "Accuracy", "Smoothness"),
                           density = c("uniform", "rw"),
                           passband = 2*pi/12,
                           optimal.minBandwidth = horizon,
                           optimal.maxBandwidth = 3*horizon){
  kernel <- match.arg(tolower(kernel)[1],
                      choices = c("biweight", "henderson", "epanechnikov", "triangular", "uniform",
                                  "triweight"))
  asymmetricCriterion = switch(tolower(asymmetricCriterion[1]),
                               timeliness = "Timeliness",
                               frequencyresponse = "FrequencyResponse",
                               accuracy = "Accuracy",
                               smoothness = "Smoothness",
                               undefined = "Undefined")
  density = match.arg(density)
  optimalBw= J("jdplus/filters/base/r/RKHSFilters")$optimalBandwidth(
    as.integer(horizon), as.integer(degree), kernel,
    asymmetricCriterion, density=="rw", passband, optimal.minBandwidth, optimal.maxBandwidth
  )
  names(optimalBw) <- sprintf("q=%i", 0:(horizon-1))
  optimalBw
}
#' Get RKHS kernel function
#' @inheritParams rkhs_filter
#' @export
rkhs_kernel <- function(kernel = c("Biweight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "Triweight"),
                        degree = 2, horizon = 6){
  kernel <- match.arg(tolower(kernel)[1],
                   choices = c("biweight", "henderson", "epanechnikov", "triangular", "uniform",
                               "triweight"))
  kernel =  switch(tolower(kernel),
    "biweight" = "BiWeight",
    "triweight" ="TriWeight",
    "uniform" = "Uniform",
    "triangular" = "Triangular",
    "epanechnikov" = "Epanechnikov",
    "henderson" = "Henderson"
  )
  kernel_fun = J("jdplus/filters/base/r/RKHSFilters")$kernel(
    kernel, as.integer(degree), as.integer(horizon)
  )$applyAsDouble

  Vectorize(function(x){
    kernel_fun(x)
  })
}
