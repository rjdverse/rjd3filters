#' Mean Square Revision Error (mmsre) filter
#'
#' Provides an asymmetric filter \eqn{[-h, p]} based on the given symmetric
#' filter.
#'
#' See Proietti, Luati, "Real time estimation in local polynomial regression
#' with application to trend-cycle analysis".
#'
#' @param sym The symmetric filter
#' @param q The horizon of the asymmetric filter (from 0 to deg(w)/2)
#' @param U Matrix U
#' @param Z Matrix Z
#' @param dz Coefficients of the linear model. The number of the
#' coefficients and the degree of the constraints define the type of the
#' linear model.
#'
#' @details
#' The asymmetric filter minimizes the mean square revision error
#' (mmsre) relative to the symmetric filter.
#' The series follows the model
#' \deqn{
#' y=U \gamma + Z \delta + \varepsilon, \quad
#' \varepsilon \sim \mathcal N(0,\sigma^2 K^{-1}).
#' }
#'
#' With \eqn{K} a set of weights (kernel), by default (`kernel = NULL`) no weight is used.
#'
#'
#' @inheritParams lp_filter
#' @examples
#' QL <- lp_filter(endpoints = "QL", ic = 3.5)
#' LC <- lp_filter(endpoints = "LC", ic = 3.5)
#' DAF <- lp_filter(endpoints = "DAF")
#' h6 <- QL[, "q=6"]
#' # To reproduce DAF filter
#' mmsreFilter(
#'   sym = h6, q = 0,
#'   U = polynomial_matrix(l = - 6, d0 = 0, d1 = 3),
#'   kernel = "Henderson"
#' )
#' DAF[, "q=0"]
#' # To reproduce QL filter
#' mmsreFilter(
#'   sym = h6, q = 1,
#'   dz = 2 / (sqrt(pi) * 3.5),
#'   U = polynomial_matrix(l = -6, d0 = 0, d1 = 1),
#'   Z = polynomial_matrix(l = -6, d0 = 2, d1 = 2)
#' )
#' QL[, "q=1"]
#'
#' # Or using the polynomial kernel
#' mmsreFilter(
#'   sym = h6, q = 2,
#'   # we multiply by the square root of the inverse of weights (1/13)
#'   # to get the same result as the QL filter
#'   dz = 2 / (sqrt(pi) * 3.5) * (sqrt(13)),
#'   U = polynomial_matrix(l = -6, d0 = 0, d1 = 0),
#'   Z = polynomial_matrix(l = -6, d0 = 1, d1 = 1)
#' )
#' LC[, "q=2"]
#' @export
mmsreFilter <- function(
    sym, q, U, Z = NULL, dz = NULL,
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
  if (is.null(dz))
    dz <- numeric()

  jkernel <- .r2jd_kernel(kernel, abs(.jcall(jsym, "I", "getLowerBound")))
  jf <- .jcall(
    "jdplus/toolkit/base/core/math/linearfilters/AsymmetricFiltersFactory",
    "Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;",
    "mmsreFilter",
    jsym,
    as.integer(q),
    .r2jd_fast_matrix(U),
    .r2jd_fast_matrix(Z),
    .jarray(dz),
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
# dz = 3.5
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
#   "mmsreFilter",
#   jsym,
#   as.integer(q),
#   .r2jd_fast_matrix(U),
#   .r2jd_fast_matrix(Z),
#   .jarray(2 / (sqrt(pi) * dz)),
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
