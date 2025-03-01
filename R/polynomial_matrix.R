#' Create polynomial matrix
#'
#' Create polynomial matrix used in local polynomial regression (see details).
#'
#' @param l,u lower bound (usually negative) and upper bound (usually positive) of the polynomial matrix.
#' @param d0,d1 lower and polynomial degree of the polynomial matrix.
#'
#' @details
#' `polynomial_matrix()` computes the following matrix
#' \deqn{
#' \begin{pmatrix}
#' (l)^{d_0} & (l)^{d_0+1} & \cdots&(l)^{d_1}\\
#' (l+1)^{d_0} & (l+1)^{d_0+1} & \cdots&(l+1)^{d_1} \\
#' \vdots & \vdots & \cdots & \vdots \\
#' (p)^{d_0} & (p)^{d_0+1} & \cdots&(p)^{d_1}
#' \end{pmatrix}
#' }
#'
#' @examples
#' # For example to reproduce DAF filters
#' daf <- lp_filter(horizon = 6, endpoints = "DAF")
#' q <- 0
#' X <- polynomial_matrix(l = -6, u = q, d0 = 0, d1 = 3)
#' K <- diag(sapply(-6:q, function(i) get_kernel(horizon = 6)[i]))
#' e_1 <- c(1, 0, 0, 0)
#' q0 <- K %*% X %*% solve(t(X) %*% K %*% X, e_1)
#' q0
#' daf[, "q=0"]
#' @export
polynomial_matrix <- function(l, u = abs(l), d0 = 0, d1 = 3) {
  .jd2r_matrix(
    .jcall(
      "jdplus/toolkit/base/core/math/linearfilters/LocalPolynomialFilters",
      "Ljdplus/toolkit/base/core/math/matrices/FastMatrix;",
      "z",
      .jnull("jdplus/toolkit/base/core/math/matrices/FastMatrix"),
      as.integer(l), as.integer(u),
      as.integer(d0), as.integer(d1))
  )
}
