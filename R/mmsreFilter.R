#' Mean Square Revision Error (mmsre) filter
#'
#' Provides an asymmetric filter based on the given reference
#' filter (usually symmetric) minimizing the mean square revision error.
#'
#' @param ref_filter The reference filter (a [moving_average()] object).
#' @param q The horizon of the asymmetric filter.
#' @param U Matrix of the constraints.
#' @param Z Matrix of the bias (can be `NULL`).
#' @param delta Coefficients of the linear model.
#' @param kernel The kernel used for weighting factors, by default, no weight is used.
#' See [lp_filter()] for the available kernels.
#'
#' @details
#' The asymmetric filter \eqn{\boldsymbol v=(v_{-h},\dots,v{q})'} minimizes the mean square revision error
#' (mmsre) relative to the reference filter \eqn{\boldsymbol \theta=(\theta_{-h},\dots,\theta_{h'})'}.
#' The series follows the model
#' \deqn{
#' \boldsymbol y=\boldsymbol U \boldsymbol \gamma +
#' \boldsymbol Z \boldsymbol \delta + \boldsymbol \varepsilon, \quad
#' \boldsymbol \varepsilon \sim \mathcal N(0,\sigma^2 \boldsymbol K^{-1}).
#' }
#'
#' With \eqn{K} a set of weights (kernel), by default (`kernel = NULL`) no weight is used.
#' The matrix \eqn{U} represents the constraints of the symmetric filter (usually polynomials preservations), \eqn{\boldsymbol \theta},
#' imposed to the asymmetric filter, \eqn{\boldsymbol v}.
#' Partitionning the matrix \eqn{\boldsymbol U=\begin{pmatrix} \boldsymbol U_p' & \boldsymbol U_f'\end{pmatrix}'}
#' with \eqn{\boldsymbol U_p} the first \eqn{h+q+1} rows and \eqn{\boldsymbol U_f} the remaining rows, the constraints are
#' \eqn{\boldsymbol U_p'\boldsymbol v=\boldsymbol U'\boldsymbol \theta}.
#'
#' The matrix \eqn{\boldsymbol Z} represents the bias of the asymmetric filter: usually constraints imposed to the symmetric filter but not to the asymmetric filter.
#'
#'
#'
#'
#' @inheritParams lp_filter
#' @examples
#' QL <- lp_filter(endpoints = "QL", ic = 3.5)
#' LC <- lp_filter(endpoints = "LC", ic = 3.5)
#' DAF <- lp_filter(endpoints = "DAF")
#' h6 <- QL[, "q=6"]
#' # To reproduce DAF filter
#' mmsre_filter(
#'   ref_filter = h6, q = 0,
#'   U = polynomial_matrix(l = - 6, d0 = 0, d1 = 3),
#'   kernel = "Henderson"
#' )
#' DAF[, "q=0"]
#' # To reproduce QL filter
#' mmsre_filter(
#'   ref_filter = h6, q = 1,
#'   delta = 2 / (sqrt(pi) * 3.5),
#'   U = polynomial_matrix(l = -6, d0 = 0, d1 = 1),
#'   Z = polynomial_matrix(l = -6, d0 = 2, d1 = 2)
#' )
#' QL[, "q=1"]
#'
#' # Or using the Uniform kernel
#' mmsre_filter(
#'   ref_filter = h6, q = 2,
#'   # we multiply by the square root of the inverse of weights (1/13)
#'   # to get the same result as the QL filter
#'   delta = 2 / (sqrt(pi) * 3.5) * (sqrt(13)),
#'   U = polynomial_matrix(l = -6, d0 = 0, d1 = 0),
#'   Z = polynomial_matrix(l = -6, d0 = 1, d1 = 1),
#'   kernel = "Uniform"
#' )
#' LC[, "q=2"]
#' @seealso [lp_filter()].
#' @references Proietti, Tommaso and Alessandra Luati (2008). “Real time estimation in local polynomial regression, with application to trend-cycle analysis”.
#' @export
mmsre_filter <- function(
    ref_filter, q, U, Z = NULL, delta = NULL,
    kernel = NULL,
    tweight = 0, passband = pi/12){
  jref <- .jcast(.ma2jd(ref_filter), "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter")
  if (is.null(delta))
    delta <- numeric()

  jkernel <- .r2jd_kernel(kernel, abs(.jcall(jref, "I", "getLowerBound")))
  jf <- .jcall(
    "jdplus/toolkit/base/core/math/linearfilters/AsymmetricFiltersFactory",
    "Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;",
    "mmsreFilter",
    jref,
    as.integer(q),
    .r2jd_fast_matrix(U),
    .r2jd_fast_matrix(Z),
    .jarray(delta),
    jkernel,
    as.numeric(tweight),
    as.numeric(passband)
  )
  return(.jd2ma(jf))
}

.r2jd_fast_matrix <- function(s){
  if (is.null(s))
    return(.jnull("jdplus/toolkit/base/core/math/matrices/FastMatrix"))

  .jcall(
    "jdplus/toolkit/base/core/math/matrices/FastMatrix",
    "Ljdplus/toolkit/base/core/math/matrices/FastMatrix;",
    "of",
    rjd3toolkit::.r2jd_matrix(s))
}

