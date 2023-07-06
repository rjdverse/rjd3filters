#' Estimation of a filter using the Fidelity-Smoothness-Timeliness criteria
#'
#' @param lags Lags of the filter (should be positive).
#' @param leads Leads of the filter (should be positive or 0).
#' @param pdegree Local polynomials preservation: max degree.
#' @param smoothness.weight Weight for the smoothness criterion (in \eqn{[0, 1]}).
#' @param smoothness.degree Degree of the smoothness criterion (3 for Henderson).
#' @param timeliness.weight Weight for the Timeliness criterion (in \eqn{[0, 1[}). \code{sweight+tweight} should be in \eqn{[0,1]}.
#' @param timeliness.passband Passband for the timeliness criterion (in radians). The phase effect is computed in \eqn{[0, passband]}.
#' @param timeliness.antiphase boolean indicating if the timeliness should be computed analytically (\code{TRUE}) or numerically (\code{FALSE}).
#'
#' @return An object of class \code{"fst_filter"}, which is a list of 5 elements:\itemize{
#' \item{\code{"internal"}}{Java object used for internal functions}
#' \item{\code{"filters.coef"}}{The coefficients of the selected filter}
#' \item{\code{"filters.gain"}}{The gain function between 0 and pi (601 observations)}
#' \item{\code{"filters.phase"}}{The phase function between 0 and pi (601 observations)}
#' \item{\code{"criteria"}}{The value of the tree criteria}
#' }
#'
#' @examples
#' filter <- fst_filter(lags = 6, leads = 0)
#' filter
#' @references Grun-Rehomme, Michel, Fabien Guggemos, and Dominique Ladiray (2018). “Asymmetric Moving Averages Minimizing Phase Shift”. In: Handbook on Seasonal Adjustment.
#' @export
fst_filter<-function(lags = 6, leads = 0, pdegree = 2,
                    smoothness.weight = 1, smoothness.degree = 3, timeliness.weight = 0,
                    timeliness.passband = pi/6, timeliness.antiphase = TRUE){
  jobj<-.jcall("jdplus/filters/base/r/FiltersToolkit", "Ljdplus/filters/base/r/FiltersToolkit$FSTResult;",
               "fstfilter", as.integer(lags), as.integer(leads), as.integer(pdegree), smoothness.weight, as.integer(smoothness.degree),
               timeliness.weight, timeliness.passband, as.logical(timeliness.antiphase))
  jfilter <- .jcall(jobj, "Ljdplus/toolkit/base/core/math/linearfilters/FiniteFilter;", "getFilter")

  return(.jd2ma(jfilter))
}

#' FST criteria
#'
#' Compute the Fidelity, Smoothness and Timeliness (FST) criteria
#'
#' @param weights either a `"moving_average"` or a numeric vector containing weights.
#' @param lags Lags of the moving average (when `weights` is not a `"moving_average"`).
#' @param passband Passband threshold for timeliness criterion.
#'
#' @return The values of the 3 criteria, the gain and phase of the associated filter.
#' @examples
#' filter <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "LC")
#' fst(filter[, "q=0"])
#' @references Grun-Rehomme, Michel, Fabien Guggemos, and Dominique Ladiray (2018). “Asymmetric Moving Averages Minimizing Phase Shift”. In: Handbook on Seasonal Adjustment.
#'
#' @export
fst<-function(weights, lags, passband=pi/6){
  if (is.moving_average(weights)) {
    lags <- lower_bound(weights)
    weights <- coef(weights)
  }
  jobj<-.jcall("jdplus/filters/base/r/FiltersToolkit", "Ljdplus/filters/base/r/FiltersToolkit$FSTResult;", "fst",
               weights, as.integer(lags), passband)
  criteria<-.jcall(jobj, "[D", "getCriterions")
  names(criteria) <- c("Fidelity", "Smoothness", "Timeliness")
  return(criteria)
}

#' Accuracy/smoothness/timeliness criteria through spectral decomposition
#'
#'
#' @param sweights `moving_average` object or weights of the symmetric filter (from 0 to n or -n to n).
#' @param aweights `moving_average` object or weights of the asymmetric filter (from -n to m).
#' @param density hypothesis on the spectral density: \code{"uniform"} (= white woise, the default) or  \code{"rw"} (= random walk).
#' @param passband passband threshold.
#'
#' @return The criteria
#' @examples
#' filter <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "LC")
#' sweights <- filter[,"q=6"]
#' aweights <- filter[,"q=0"]
#' mse(sweights, aweights)
#' @references Wildi, Marc and McElroy, Tucker (2019). “The trilemma between accuracy, timeliness and smoothness in real-time signal extraction”. In: International Journal of Forecasting 35.3, pp. 1072–1084.
#' @export
mse<-function(sweights, aweights, density=c("uniform", "rw"), passband = pi/6){
  if (is.moving_average(aweights))
    aweights <- coef(aweights)

  if (is.moving_average(sweights)) {
    if (lower_bound(sweights) < 0) {
      # we asume sweights were specify from [-n to n] instead of [0,n]
      sweights <- coef(sweights)[seq(lower_bound(sweights), -1)]
    } else {
      sweights <- coef(sweights)
    }
  } else if(length(sweights)>length(aweights)){
    # we asume sweights were specify from [-n to n] instead of [0,n]
    n <- (length(sweights)-1)/2
    sweights <- sweights[-seq_len(n)]
  }
  spectral = match.arg(density)
  rslt<-.jcall("jdplus/filters/base/r/FiltersToolkit", "[D", "mseDecomposition",
               sweights, aweights, spectral, passband)
  return (c(accuracy=rslt[1], smoothness=rslt[2], timeliness=rslt[3], residual=rslt[4]))
}
