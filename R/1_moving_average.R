#'@importFrom methods is
setClass("moving_average",
         slots = c(coefficients = "numeric",
                   lower_bound = "numeric",
                   upper_bound = "numeric")
)

#' Operations on Filters
#'
#' Manipulation of [moving_average()] or [finite_filters()] objects
#'
#' @param x,e1,e2 object
#' @param i,j,value indices specifying elements to extract or replace and the new value
#'
#' @param ...,drop,na.rm other parameters.
#' @name filters_operations
NULL

#' Manipulation of moving averages
#'
#' @param x vector of coefficients.
#' @param lags integer indicating the number of lags of the moving average.
#' @param trailing_zero,leading_zero boolean indicating wheter to remove leading/trailing zero and NA.
#' @param s seasonal period for the \code{to_seasonal()} function.
#' @param object `moving_average` object.
#'
#' @details
#' A moving average is defined by a set of coefficient \eqn{\boldsymbol \theta=(\theta_{-p},\dots,\theta_{f})'}
#' such all time series \eqn{X_t} are transformed as:
#' \deqn{
#' M_{\boldsymbol\theta}(X_t)=\sum_{k=-p}^{+f}\theta_kX_{t+k}=\left(\sum_{k=-p}^{+f}\theta_kB^{-k}\right)X_{t}
#' }
#' The integer \eqn{p} is defined by the parameter \code{lags}.
#'
#' The function `to_seasonal()` transforms the moving average \eqn{\boldsymbol \theta} to:
#' \deqn{
#' M_{\boldsymbol\theta'}(X_t)=\sum_{k=-p}^{+f}\theta_kX_{t+ks}=\left(\sum_{k=-p}^{+f}\theta_kB^{-ks}\right)X_{t}
#' }
#'
#' @examples
#' y <- retailsa$AllOtherGenMerchandiseStores
#' e1 <- moving_average(rep(1,12), lags = -6)
#' e1 <- e1/sum(e1)
#' e2 <- moving_average(rep(1/12, 12), lags = -5)
#' M2X12 <- (e1 + e2)/2
#' coef(M2X12)
#' M3 <- moving_average(rep(1/3, 3), lags = -1)
#' M3X3 <- M3 * M3
#' # M3X3 moving average applied to each month
#' M3X3
#' M3X3_seasonal <- to_seasonal(M3X3, 12)
#' # M3X3_seasonal moving average applied to the global series
#' M3X3_seasonal
#'
#' def.par <- par(no.readonly = TRUE)
#' par(mai = c(0.5, 0.8, 0.3, 0))
#' layout(matrix(c(1,2), nrow = 1))
#' plot_gain(M3X3, main = "M3X3 applied to each month")
#' plot_gain(M3X3_seasonal, main = "M3X3 applied to the global series")
#' par(def.par)
#'
#' # To apply the moving average
#' t <- y * M2X12
#' # Or use the filter() function:
#' t <- filter(y, M2X12)
#' si <- y - t
#' s <- si * M3X3_seasonal
#' # or equivalently:
#' s_mm <- M3X3_seasonal * (1 - M2X12)
#' s <- y * s_mm
#' plot(s)
#' @export
moving_average <- function(x, lags = -length(x), trailing_zero = FALSE, leading_zero = FALSE){
  if (inherits(x, "moving_average"))
    return (x)
  x <- as.numeric(x)
  if (trailing_zero)
    x <- rm_trailing_zero_or_na(x)
  if (leading_zero){
    new_x <- rm_leading_zero_or_na(x)
    lags <- lags - (length(new_x) - length(x))
    x <- new_x
  }
  upper_bound <- lags + length(x) -1
  # remove 1 if it is >= 0 (central term)
  # upper_bound = upper_bound - (upper_bound >= 0)

  names(x) <- coefficients_names(lags,
                                 upper_bound)
  res <- new("moving_average",
             coefficients = x, lower_bound = lags,
             upper_bound = upper_bound)
  res
}
.jd2ma <- function(jobj, trailing_zero = FALSE){
  x <- .jcall(jobj, "[D", "weightsToArray")
  lags <- .jcall(jobj, "I", "getLowerBound")
  moving_average(x, lags, trailing_zero = trailing_zero)
}
.ma2jd <- function(x){
  lags <- lower_bound(x)
  coefs <- as.numeric(coef(x))
  if (length(x) == 1){
    coefs <- .jarray(coefs)
  }
  .jcall("jdplus/toolkit/base/core/math/linearfilters/FiniteFilter",
         "Ljdplus/toolkit/base/core/math/linearfilters/FiniteFilter;",
         "of", coefs,
         as.integer(lags))
}
#' @rdname moving_average
#' @export
is.moving_average <- function(x){
  is(x, "moving_average")
}
#' @importFrom stats coef coefficients
#' @export
coef.moving_average <- function(object, ...){
  coefs <- object@coefficients
  return(coefs)
}
#' @rdname moving_average
#' @export
is_symmetric <- function(x){
  # .jcall(.ma2jd(x), "Z", "isSymmetric")
  (upper_bound(x) == (-lower_bound(x))) &&
    isTRUE(all.equal(coef(x), rev(coef(x)), check.attributes = FALSE))
}
#' @rdname moving_average
#' @export
upper_bound <- function(x){
  x@upper_bound
}
#' @rdname moving_average
#' @export
lower_bound <- function(x){
  x@lower_bound
}
#' @rdname moving_average
#' @export
mirror <- function(x){
  .jd2ma(.jcall(.ma2jd(x), "Ljdplus/toolkit/base/core/math/linearfilters/FiniteFilter;", "mirror"))
}
#' @method rev moving_average
#' @rdname moving_average
#' @export
rev.moving_average <- function(x){
  mirror(x)
}
#' @rdname moving_average
#' @export
length.moving_average <- function(x){
  length(coef(x))
}
#' @rdname moving_average
#' @export
to_seasonal <- function(x, s){
  UseMethod("to_seasonal", x)
}
#' @export
to_seasonal.default <- function(x, s){
  lb <- lower_bound(x)
  up <- upper_bound(x)
  coefs <- coef(x)
  new_coefs <- c(unlist(lapply(coefs[-length(x)],
                               function(x){
                                 c(x, rep(0, s - 1))
                               })),
                 coefs[length(x)])
  moving_average(new_coefs, lb * s)
}

#' @rdname filters_operations
#' @export
sum.moving_average <- function(..., na.rm = FALSE){
  sum(
    unlist(lapply(list(...),
                  function(x) sum(coef(x),na.rm = na.rm)
    )
    )
  )
}
#' @rdname filters_operations
#' @export
setMethod("[",
          signature(x = "moving_average",
                    i = "numeric"),
          function(x, i) {
            coefs <- coef(x)
            indices <- seq_along(coefs)[i]
            coefs[-indices] <- 0
            if (all(coefs == 0))
              return (moving_average(0, lags = lower_bound(x) + indices - 1))

            moving_average(coefs, lags = lower_bound(x),
                           leading_zero = TRUE, trailing_zero = TRUE)
          })
#' @rdname filters_operations
#' @export
setMethod("[",
          signature(x = "moving_average",
                    i = "logical"),
          function(x, i) {
            coefs <- coef(x)
            indices <- seq_along(coefs)[i]
            coefs[!indices] <- 0
            moving_average(coefs, lags = lower_bound(x),
                           leading_zero = TRUE, trailing_zero = TRUE)
          })
#' @rdname filters_operations
#' @export
setReplaceMethod("[",
                 signature(x = "moving_average",
                           i = "ANY",
                           j = "missing",
                           value = "numeric"),
                 function(x, i, value) {
                   x@coefficients[i] <- value
                   x
                 })
#' @rdname filters_operations
#' @export
cbind.moving_average <- function(..., zero_as_na = FALSE){
  all_mm <- list(...)
  new_lb <- min(sapply(all_mm, lower_bound))
  new_ub <- max(sapply(all_mm, upper_bound))
  nb_uterms <- max(sapply(all_mm, function(x) lower_bound(x) + length(x)))
  if (zero_as_na) {
      blank_value <- NA
  } else {
      blank_value <- 0
  }
  new_mm <- lapply(all_mm, function(x){
    c(rep(blank_value, abs(new_lb - lower_bound(x))),
      coef(x),
      rep(blank_value, abs(nb_uterms - (lower_bound(x) + length(x))))
    )
  })
  new_mm <- do.call(cbind, new_mm)
  rownames(new_mm) <- coefficients_names(new_lb, new_ub)
  new_mm
}
#' @rdname filters_operations
#' @export
rbind.moving_average <- function(...){
  t(cbind(...))
}
#' @rdname filters_operations
#' @export
setMethod("+",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.toolkit.base.core.math.linearfilters.FiniteFilter")
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/toolkit/base/core/math/linearfilters/FiniteFilter;",
                           "add",
                           .jcast(.ma2jd(e1), "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter"),
                           .jcast(.ma2jd(e2), "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter"))

            .jd2ma(jobj)
          })
#' @rdname filters_operations
#' @export
setMethod("+",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 + moving_average(e2,0)
          })
#' @rdname filters_operations
#' @export
setMethod("+",
          signature(e1 = "numeric",
                    e2 = "moving_average"),
          function(e1, e2) {
            e2 + e1
          })
#' @rdname filters_operations
#' @export
setMethod("+", signature(e1 = "moving_average", e2 = "missing"), function(e1,e2) e1)
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "missing"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.toolkit.base.core.math.linearfilters.FiniteFilter")
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/toolkit/base/core/math/linearfilters/FiniteFilter;",
                           "negate",
                           .jcast(.ma2jd(e1), "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter"))
            .jd2ma(jobj)
          })
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.toolkit.base.core.math.linearfilters.FiniteFilter")
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/toolkit/base/core/math/linearfilters/FiniteFilter;",
                           "subtract",
                           .jcast(.ma2jd(e1), "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter"),
                           .jcast(.ma2jd(e2), "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter"))
            .jd2ma(jobj)
          })
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 + (- e2)
          })
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "numeric",
                    e2 = "moving_average"),
          function(e1, e2) {
            e1 + (- e2)
          })
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.toolkit.base.core.math.linearfilters.FiniteFilter")
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/toolkit/base/core/math/linearfilters/FiniteFilter;",
                           "multiply",
                           .jcast(.ma2jd(e1), "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter"),
                           .jcast(.ma2jd(e2), "jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter"))
            .jd2ma(jobj)
          })
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            if (length(e2) == 1) {
              e1 * moving_average(e2,0)
            } else {
              filter(e2, e1)
            }
          })

#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "numeric",
                    e2 = "moving_average"),
          function(e1, e2) {
            if (length(e1) == 1) {
              moving_average(e1,0) * e2
            } else {
              filter(e1, e2)
            }
          })
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e2 = "moving_average"),
          function(e1, e2) {
            filter(e1,e2)
          })
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "moving_average"),
          function(e1, e2) {
            filter(e2, e1)
          })
#' @rdname filters_operations
#' @export
setMethod("/",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 * moving_average(1/e2,0)
          })
#' @rdname filters_operations
#' @export
setMethod("^",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            if (e2 == 0) {
              moving_average(1, 0)
            } else {
              Reduce(`*`, rep(list(e1), e2))
            }
          })
#' Simple Moving Average
#'
#' A simple moving average is a moving average whose coefficients are all equal and whose sum is 1
#'
#' @param order number of terms of the moving_average
#' @inheritParams moving_average
#'
#' @examples
#' # The M2X12 moving average is computed as
#' (simple_ma(12, -6) + simple_ma(12, -5)) / 2
#' # The M3X3 moving average is computed as
#' simple_ma(3, -1) ^ 2
#' # The M3X5 moving average is computed as
#' simple_ma(3, -1) * simple_ma(5, -2)
#' @export
simple_ma <- function(order, lags = - trunc((order-1)/2)) {
  moving_average(rep(1, order), lags = lags) / order
}
#'@export
as.list.moving_average <- function(x, ...) {
  lapply(seq_along(x), function(i) x[i])
}
