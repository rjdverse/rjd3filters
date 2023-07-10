#' Direct Filter Approach
#'
#' @inheritParams localpolynomials
#' @inheritParams mse
#' @param targetfilter the weights of the symmetric target filters (by default the Henderson filter).
#' @param accuracy.weight,smoothness.weight,timeliness.weight the weight used for the
#' optimisation. The weight associated to the residual is derived so that the sum of
#' the four weights are equal to 1.
#' @export
#' @examples
#' dfa_filter(horizon = 6, degree = 0)
#' dfa_filter(horizon = 6, degree = 2)
dfa_filter <- function(horizon = 6, degree = 0,
                       density = c("uniform", "rw"),
                       targetfilter = lp_filter(horizon = horizon)[,1],
                       passband = 2*pi/12,
                       accuracy.weight = 1/3,
                       smoothness.weight = 1/3,
                       timeliness.weight = 1/3){
  density = match.arg(density)
  if (length(targetfilter) != 2*horizon + 1)
    stop("The symmetric targetfilter must be of length 2*horizon+1")
  if (is.moving_average(targetfilter)) {
    if (lower_bound(targetfilter) < 0) {
      # we asume targetfilter were specify from [-n to n] instead of [0,n]
      targetfilter <- coef(targetfilter)[seq(lower_bound(targetfilter), -1)]
    } else {
      targetfilter <- coef(targetfilter)
    }
  }
  dfa_filter = J("jdplus/filters/base/r/DFAFilters")$filterProperties(
    targetfilter,
    as.integer(horizon), as.integer(degree), density=="rw",
    passband,
    accuracy.weight, smoothness.weight, timeliness.weight
  )
  return(.jd2r_finitefilters(dfa_filter))
}
#

