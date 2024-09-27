#' Mean Square Revision Error (mmsre) filter
#'
#' Provides an asymmetric filter based on the given symmetric
#' filter minimizing the mean square revision error.
#'
#' @param sym The symmetric filter.
#' @param q The horizon of the asymmetric filter.
#' @param U Matrix of the constraints.
#' @param Z Matrix of the bias (can be `NULL`).
#' @param delta Coefficients of the linear model.
#' @param kernel The kernel used for weighting factors, by default, no weight is used.
#' See [lp_filter()] for the available kernels.
#'
#' @details
#' The asymmetric filter \eqn{\boldsymbol v=(v_{-h},\dots,v{q})'} minimizes the mean square revision error
#' (mmsre) relative to the symmetric filter \eqn{\boldsymbol \theta=(\theta_{-h},\dots,\theta_{h})'}.
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
#'   sym = h6, q = 0,
#'   U = polynomial_matrix(l = - 6, d0 = 0, d1 = 3),
#'   kernel = "Henderson"
#' )
#' DAF[, "q=0"]
#' # To reproduce QL filter
#' mmsre_filter(
#'   sym = h6, q = 1,
#'   delta = 2 / (sqrt(pi) * 3.5),
#'   U = polynomial_matrix(l = -6, d0 = 0, d1 = 1),
#'   Z = polynomial_matrix(l = -6, d0 = 2, d1 = 2)
#' )
#' QL[, "q=1"]
#'
#' # Or using the Uniform kernel
#' mmsre_filter(
#'   sym = h6, q = 2,
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
    sym, q, U, Z = NULL, delta = NULL,
    kernel = NULL,
    tweight = 0, passband = pi/12){
  if (is.moving_average(sym)) {
    sym <- coef(sym)
  } else if (is.finite_filters(sym)) {
    sym <- coef(sym@sfilter)
  }
  jsym <- .jcall(
    "jdplus/toolkit/base/core/math/linearfilters/SymmetricFilter",
    "Ljdplus/toolkit/base/core/math/linearfilters/SymmetricFilter;",
    "of", .r2jd_doubleseq(sym)
  )
  if (is.null(delta))
    delta <- numeric()

  jkernel <- .r2jd_kernel(kernel, abs(.jcall(jsym, "I", "getLowerBound")))
  jf <- .jcall(
    "jdplus/toolkit/base/core/math/linearfilters/AsymmetricFiltersFactory",
    "Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;",
    "mmsreFilter",
    jsym,
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
# library(rJava)
# devtools::load_all()
# sym <- lp_filter()[,"q=6"]
# U = polynomial_matrix(l = -6, d0 = 0, d1 = 1)
# Z = polynomial_matrix(l = -6, d0 = 2, d1 = 2)
# q=0
#
# delta = 3.5
# tweight = 0
# passband = pi/12
# kernel <- "uniform"
# jsym <- .jcall(
#   "jdplus/toolkit/base/core/math/linearfilters/SymmetricFilter",
#   "Ljdplus/toolkit/base/core/math/linearfilters/SymmetricFilter;",
#   "of", rjd3filters:::.r2jd_doubleseq(coef(sym))
# )
# jkernel <- .r2jd_kernel(kernel, abs(.jcall(jsym, "I", "getLowerBound")))
# .jcall(
#   "jdplus/toolkit/base/core/math/linearfilters/AsymmetricFiltersFactory",
#   "Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;",
#   "mmsre_filter",
#   jsym,
#   as.integer(q),
#   .r2jd_fast_matrix(U),
#   .r2jd_fast_matrix(Z),
#   .jarray(2 / (sqrt(pi) * delta)),
#   jkernel,
#   as.numeric(tweight),
#   as.numeric(passband)
# )
# lp_filter(ic = 3.5)[,"q=0"]
#
# jU = .r2jd_fast_matrix(U)
# jZ = .r2jd_fast_matrix(Z)
# w = jsym$weightsToArray()
# h = as.integer(length(w) / 2)
# nv = as.integer(h + q + 1)
# ncolu = as.integer(jU$getColumnsCount())
# DataBlock <- J("jdplus.toolkit.base.core.data.DataBlock")
# FastMatrix <- J("jdplus/toolkit/base/core/math/matrices/FastMatrix")
# wp = DataBlock$of(w, 0L, as.integer(nv));
# wf = DataBlock$of(w, nv, length(w));
# Up = jU$extract(0L, nv, 0L, ncolu);
# Uf = jU$extract(nv, jU$getRowsCount()-nv, 0L, ncolu);
# Zp = jZ$extract(0L, nv, 0L, jZ$getColumnsCount());
# Zf = jZ$extract(nv, jZ$getRowsCount() - nv, 0L, jZ$getColumnsCount());
# Q = FastMatrix$square(nv + ncolu);
# D = Q$extract(0L, nv, 0L, nv);
# D$diagonal()$set(1);
# D$diagonal()$set(1L,1);
# .jcall(Q$extract(nv, ncolu, 0L, nv), "V", "copyTranspose", Up)
# .jcall(Q$extract(0L,nv, nv, ncolu), "V", "copy", Up)
#
# Q.extract(nv, u + 1, 0, nv).copyTranspose(Up);
# Q.extract(0, nv, nv, u + 1).copy(Up);
# LocalPolynomialFilters <- J("jdplus.toolkit.base.core.math.linearfilters.LocalPolynomialFilters")
# Z = LocalPolynomialFilters$createZ(6L, 1L);
# Up = LocalPolynomialFilters$z(
#   LocalPolynomialFilters$createZ(6L, 1L),
#   -6L, 0L, 0L, 0L);
# Uf = LocalPolynomialFilters$z(
#   LocalPolynomialFilters$createZ(6L, 1L),
#   as.integer(q + 1), as.integer(h), 0L, 0L);
# Q = FastMatrix$square(nv + ncolu);
# D = Q$extract(0L, nv, 0L, nv);
# D$diagonal()$set(1);
# .jcall(Q$extract(nv, ncolu, 0L, nv), "V", "copyTranspose", Up)
# .jcall(Q$extract(0L,nv, nv, ncolu), "V", "copy", Up)
# DataBlock a = DataBlock.make(Q.getRowsCount());
# a.extract(nv, ncolu + 1).product(wf, Uf.columnsIterator());
# .jmethods("jdplus/toolkit/base/core/math/linearfilters/SymmetricFilter","of")
#
# Q.extract(nv, u + 1, 0, nv).copyTranspose(Up);
# Q.extract(0, nv, nv, u + 1).copy(Up);
