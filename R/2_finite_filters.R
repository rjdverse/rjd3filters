setClass("finite_filters",
         slots = c(sfilter = "moving_average",
                   lfilters = "list",
                   rfilters = "list")
)
#' Manipulating Finite Filters
#'
#' @param sfilter the symmetric filter ([moving_average()] object) or
#'  a `matrix` or `list` with all the coefficients.
#' @param rfilters the right filters (used on the last points).
#' @param lfilters the left filters (used on the first points).
#' @param first_to_last boolean indicating if the first element of `rfilters` is the
#' first asymmetric filter (when only one observation is missing) or the last one (real-time estimates).
#' @param object `finite_filters` object.
#' @param x object to test the class.
#'
#' @examples
#' ff_lp <- lp_filter()
#' ff_simple_ma <- finite_filters(moving_average(c(1, 1, 1), lags = -1)/3,
#'                rfilters = list(moving_average(c(1, 1), lags = -1)/2))
#' ff_lp
#' ff_simple_ma
#' ff_lp * ff_simple_ma
#'
#' @export
finite_filters <- function(sfilter,
                           rfilters = NULL,
                           lfilters = NULL,
                           first_to_last = FALSE){
  UseMethod("finite_filters", sfilter)
}
#' @export
finite_filters.moving_average <- function(sfilter,
                                          rfilters = NULL,
                                          lfilters = NULL,
                                          first_to_last = FALSE){
  if (is.null(lfilters) & !is.null(rfilters)) {
    if (first_to_last) {
      rfilters <- rev(rfilters)
    }
    lfilters <- rev(lapply(rfilters, rev.moving_average))
  } else if (!is.null(lfilters) & is.null(rfilters)) {
    if (!first_to_last) {
      lfilters <- rev(lfilters)
    }
    rfilters <- rev(lapply(lfilters, rev.moving_average))
  } else if (is.null(lfilters) & is.null(rfilters)) {
    rfilters <- lfilters <- list()

  }
  res <- new("finite_filters",
             sfilter = sfilter, lfilters = lfilters,
             rfilters = rfilters)
  res
}

#' @export
finite_filters.list <- function(sfilter,
                                rfilters = NULL,
                                lfilters = NULL,
                                first_to_last = FALSE){
  lags <- length(sfilter)-1

  all_f <- lapply(sfilter,
                   function(x){
                     moving_average(rm_trailing_zero_or_na(x), -lags)
                   }
  )
  if (first_to_last)
    all_f <- rev(all_f)
  sfilter <- all_f[[1]]
  rfilters <- all_f[-1]
  finite_filters(sfilter = sfilter, rfilters = rfilters)
}

#' @export
finite_filters.matrix <- function(sfilter,
                                  rfilters = NULL,
                                  lfilters = NULL,
                                  first_to_last = FALSE){
  coefs <- lapply(1:ncol(sfilter), function(i) sfilter[,i])
  finite_filters(coefs, first_to_last = first_to_last)
}

#' @rdname finite_filters
#' @export
is.finite_filters <- function(x){
    is(x, "finite_filters")
}

#' @export
.jd2r_finitefilters <- function(jf, first_to_last){
    jf<-.jcast(jf, "jdplus.toolkit.base.core.math.linearfilters/IFiltering")
    if (! is.jnull(jf)) {
      jsfilter <- .jcall(jf, "Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;", "centralFilter")
      jrfilter <- .jcall(jf, "[Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;", "rightEndPointsFilters")
      jlfilter <- .jcall(jf, "[Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;", "leftEndPointsFilters")

      sfilter <- .jd2ma(jsfilter)
      rfilters <- lapply(jrfilter, .jd2ma)
      lfilters <- rev(lapply(jlfilter, .jd2ma))

      if (missing(first_to_last)) {
        if (all(diff(sapply(lfilters, length)) <= 0)) {
          lfilters <- rev(lfilters)
          rfilters <- rev(rfilters)
        }
      } else {
        if (first_to_last) {
          lfilters <- rev(lfilters)
          rfilters <- rev(rfilters)
        }
      }

      finite_filters(sfilter = sfilter,
                     rfilters = rfilters,
                     lfilters = lfilters)
    } else {
      NULL
    }
#  if (.jinstanceof(jf, "jdplus/x12plus/base/core/X11SeasonalFiltersFactory$AnyFilter")) {
#    jsfilter <- .jcall(jf, "Ljdplus/toolkit/base/core/math/linearfilters/SymmetricFilter;", "symmetricFilter")
#    jlfilter <- .jcall(jf, "[Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;", "leftEndPointsFilters")
#    jrfilter <- .jcall(jf, "[Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;", "rightEndPointsFilters")

#    sfilter = .jd2ma(jsfilter)
#    rfilters = lapply(jrfilter, .jd2ma)
#    lfilters = rev(lapply(jlfilter, .jd2ma))
#  } else if (.jinstanceof(jf, "jdplus/toolkit/base/core/math/linearfilters/FiltersToolkit$FiniteFilters")) {
#    jsfilter <- .jcall(jf, "Ljdplus/toolkit/base/core/math/linearfilters/SymmetricFilter;", "getFilter")
#    jrfilter <- .jcall(jf, "[Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;", "getAfilters")
#    if (!first_to_last) # lp_filter
#      jrfilter <- rev(jrfilter)

#    while (is.jnull(jrfilter[[length(jrfilter)]])) { # DFA
#      jrfilter <- jrfilter[-length(jrfilter)]
#    }
#    sfilter <- .jd2ma(jsfilter)
#    rfilters <- lapply(jrfilter, .jd2ma)
#    lfilters <- NULL
#  } else if (.jinstanceof(jf, "jdplus/filters/base/core/filters/Filtering")) {
#    jsfilter <- .jcall(jf, "Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;", "centralFilter")
#    jrfilter <- .jcall(jf, "[Ljdplus/toolkit/base/core/math/linearfilters/IFiniteFilter;", "rightEndPointsFilters")

#    sfilter <- .jd2ma(jsfilter)
#    rfilters <- lapply(jrfilter, .jd2ma)
#    lfilters <- NULL
#  }

}
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "moving_average"),
          function(e1, e2) {
            sym <- e1@sfilter * e2
            rfilters <- lapply(e1@rfilters, `*`, e2)
            lfilters <- lapply(e1@lfilters, `*`, e2)

            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)
          })
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "finite_filters"),
          function(e1, e2) {

            new_ub <- e1@upper_bound + e2@sfilter@upper_bound
            new_lb_sym <- e1@lower_bound + e2@sfilter@lower_bound
            new_lb <- e1@lower_bound + -length(e2@rfilters)

            new_e2 <- c(e2@lfilters,
                        rep(list(e2@sfilter), length(e1)),
                        e2@rfilters)
            new_e1 <- rep(list(e1), length(new_e2))
            new_e2 <- lapply(1:length(new_e2), function(i){
              new_e2[[i]] * moving_average(1,
                                           lags = (new_lb + (i - 1) + e1@lower_bound * (new_lb == new_lb_sym)))
            })
            new_e1 <- lapply(1:length(new_e1), function(i){
              new_e1[[i]] * moving_average(1, lags = (new_lb + (i - 1)))
            })
            all_f <- t(do.call(cbind,c(new_e1, new_e2)))
            mat_e1 <- all_f[seq_along(new_e1),]
            mat_e2 <- all_f[-seq_along(new_e1),]
            new_mat <- (mat_e1[, seq_along(new_e2)] %*% mat_e2)[seq_len(1 + length(e2@lfilters) + length(e2@rfilters)),]

            max_lags <- min(sapply(new_e1, lower_bound), sapply(new_e2, lower_bound))

            # i_to_remove = seq_len(-(max_lags - new_lb))
            sym_mat <- new_mat[(nrow(new_mat)+1)/2, ]
            sym <- moving_average(sym_mat,
                                  lags = new_lb_sym, leading_zero = FALSE, trailing_zero = TRUE)
            rfilters <- new_mat[-(1:((nrow(new_mat)+1)/2)),, drop = FALSE]
            rfilters <- lapply(1:nrow(rfilters),function(i){
              moving_average(rfilters[i,-seq_len(i)],
                             lags = new_lb_sym, leading_zero = FALSE, trailing_zero = TRUE)
            })

            lfilters <- new_mat[(1:((nrow(new_mat)-1)/2)), , drop = FALSE]
            lfilters <- lapply(1:nrow(lfilters),function(i){
              moving_average(lfilters[i,],
                             lags = new_lb_sym + (nrow(lfilters) - i) + 1, leading_zero = FALSE, trailing_zero = TRUE)#why -1 ?
            })
            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)

          })
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
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
          signature(e1 = "ANY",
                    e2 = "finite_filters"),
          function(e1, e2) {
            filter(e1, e2)
          })
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "ANY"),
          function(e1, e2) {
            e2 * e1
          })
#' @rdname filters_operations
#' @export
setMethod("+",
          signature(e1 = "numeric",
                    e2 = "finite_filters"),
          function(e1, e2) {
            e2 + e1
          })
#' @rdname filters_operations
#' @export
setMethod("+",
          signature(e1 = "finite_filters",
                    e2 = "moving_average"),
          function(e1, e2) {
            e1@sfilter <- e1@sfilter + e2
            e1@lfilters <- lapply(e1@lfilters, `+`, e2)
            e1@rfilters <- lapply(e1@rfilters, `+`, e2)
            e1
          })
#' @rdname filters_operations
#' @export
setMethod("+",
          signature(e1 = "moving_average",
                    e2 = "finite_filters"),
          function(e1, e2) {
            e2 + e1
          })
#' @rdname filters_operations
#' @export
setMethod("+", signature(e1 = "finite_filters", e2 = "missing"), function(e1,e2) e1)
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "finite_filters",
                    e2 = "missing"),
          function(e1, e2) {
            e1@sfilter <- - e1@sfilter
            e1@lfilters <- lapply(e1@lfilters, `-`)
            e1@rfilters <- lapply(e1@rfilters, `-`)
            e1
          })
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "finite_filters",
                    e2 = "moving_average"),
          function(e1, e2) {
            e1@sfilter <- e1@sfilter - e2
            e1@lfilters <- lapply(e1@lfilters, `-`, e2)
            e1@rfilters <- lapply(e1@rfilters, `-`, e2)
            e1
          })
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "finite_filters"),
          function(e1, e2) {
            e1 + (- e2)
          })
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 - moving_average(e2,0)
          })
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "numeric",
                    e2 = "finite_filters"),
          function(e1, e2) {
            moving_average(e1,0) - e2
          })

#' @rdname filters_operations
#' @export
setMethod("/",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            e1@sfilter <- e1@sfilter / e2
            e1@lfilters <- lapply(e1@lfilters, `/`, e2)
            e1@rfilters <- lapply(e1@rfilters, `/`, e2)
            e1
          })
#' @method as.matrix finite_filters
#' @export
as.matrix.finite_filters <- function(x, sfilter = TRUE, rfilters = TRUE, lfilters = FALSE, zero_as_na = FALSE, ...) {
  sfilter_s <- rfilters_s <- lfilters_s <-
    index_s <- index_r <- index_l <- NULL
  if (!any(sfilter, rfilters, lfilters))
    return(NULL)
  if (sfilter) {
    sfilter_s <- list(x@sfilter)
    index_s <- length(x@rfilters)
  }
  if (lfilters && length(x@lfilters) > 0) {
    lfilters_s <- x@lfilters
    index_l <- seq(0, -(length(x@lfilters) - 1))
  }
  if (rfilters && length(x@rfilters) > 0) {
    rfilters_s <- x@rfilters
    index_r <- seq(length(x@rfilters) - 1, 0)
  }
  mat <- do.call(cbind, c(lfilters_s, sfilter_s, rfilters_s,
                          list(zero_as_na = zero_as_na)))
  colnames(mat) <- sprintf("q=%i", c(index_l, index_s, index_r))
  mat
}
#' @rdname filters_operations
#' @export
setMethod("^",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            Reduce(`*`, rep(list(e1), e2))
          })
#' @rdname filters_operations
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "finite_filters"),
          function(e1, e2) {
            new_ub <- length(e1@rfilters) + length(e2@rfilters)
            new_lb <- length(e1@lfilters) + length(e2@lfilters)
            new_lb_sym <- e1@sfilter@lower_bound + e2@sfilter@lower_bound

            new_e1 <- c(e1@lfilters,
                        rep(list(e1@sfilter), 1 + (new_ub - length(e1@rfilters)) + (new_lb - length(e1@lfilters))),
                        e1@rfilters)
            new_e2 <- c(e2@lfilters,
                        rep(list(e2@sfilter), 1 + (new_ub - length(e2@rfilters)) + (new_lb - length(e2@lfilters))),
                        e2@rfilters)

            new_e1 <- lapply(1:length(new_e1), function(i){
              new_e1[[i]] * moving_average(1, lags = (-new_lb + (i - 1)))
            })
            new_e2 <- lapply(1:length(new_e2), function(i){
              new_e2[[i]] * moving_average(1, lags = (-new_lb  + (i - 1)))
            })
            all_f <- t(do.call(cbind,c(new_e1, new_e2)))

            mat_e1 <- all_f[seq_along(new_e1),]
            mat_e2 <- all_f[-seq_along(new_e1),]
            new_mat <- (mat_e1[, seq_along(new_e2)] %*% mat_e2)

            max_lags <- min(sapply(new_e1, lower_bound), sapply(new_e2, lower_bound))

            sym_mat <- new_mat[(nrow(new_mat)+1)/2,]
            sym <- moving_average(sym_mat,
                                  lags = new_lb_sym, leading_zero = TRUE, trailing_zero = TRUE)
            rfilters <- new_mat[-(1:((nrow(new_mat)+1)/2)), , drop = FALSE]
            rfilters <- lapply(1:nrow(rfilters),function(i){
              moving_average(rfilters[i,],
                             lags = new_lb_sym - i, leading_zero = TRUE, trailing_zero = TRUE)
            })
            #
            lfilters <- new_mat[(1:((nrow(new_mat)-1)/2)), , drop = FALSE]
            lfilters <- lapply(1:nrow(lfilters),function(i){
              moving_average(lfilters[i,],
                             lags = new_lb_sym + (nrow(lfilters) - i) + 1, leading_zero = TRUE, trailing_zero = TRUE)#why -1 ?
            })
            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)
          })
#' @rdname filters_operations
#' @export
setMethod("+",
          signature(e1 = "finite_filters",
                    e2 = "finite_filters"),
          function(e1, e2) {

            sfilter <- e1@sfilter + e2@sfilter

            n_rfilter <- upper_bound(e1@sfilter) + upper_bound(e2@sfilter)
            n_lfilter <- lower_bound(e1@sfilter) + lower_bound(e2@sfilter)
            n_rfilter <- max(n_rfilter, 0)
            n_lfilter <- abs(min(n_lfilter, 0))

            e1_lfilters <- c(e1@lfilters,
                             rep(list(e1@sfilter),
                                 max(-(lower_bound(e2@sfilter) - lower_bound(e1@sfilter)), 0))
            )
            e2_lfilters <- c(e2@lfilters,
                             rep(list(e2@sfilter),
                                 max(-(lower_bound(e1@sfilter) - lower_bound(e2@sfilter)), 0))
            )
            e1_rfilters <- c(rep(list(e1@sfilter),
                                 max(upper_bound(e2@sfilter) - upper_bound(e1@sfilter), 0)),
                             e1@rfilters
            )
            e2_rfilters <- c(rep(list(e2@sfilter),
                                 max(upper_bound(e1@sfilter) - upper_bound(e2@sfilter), 0)),
                             e2@rfilters
            )
            e1_lfilters_f <- c(e1_lfilters, rep(list(0),
                                                max(length(e2_lfilters) - length(e1_lfilters), 0)))
            e2_lfilters_f <- c(e2_lfilters, rep(list(0),
                                                max(length(e1_lfilters) - length(e2_lfilters), 0)))
            e1_rfilters_f <- c(e1_rfilters, rep(list(0),
                                                max(length(e2_rfilters) - length(e1_rfilters), 0)))
            e2_rfilters_f <- c(e2_rfilters, rep(list(0),
                                                max(length(e1_rfilters) - length(e2_rfilters), 0)))

            lfilters <- mapply(`+`, e1_lfilters_f, e2_lfilters_f)
            rfilters <- mapply(`+`, e1_rfilters_f, e2_rfilters_f)
            finite_filters(sfilter = sfilter, rfilters = rfilters, lfilters = lfilters)
          })
#' @rdname filters_operations
#' @export
setMethod("-",
          signature(e1 = "finite_filters",
                    e2 = "finite_filters"),
          function(e1, e2) {
            e1 + (-e2)
          })
#' @rdname filters_operations
#' @export
setMethod("[",
          signature(x = "finite_filters",
                    i = "missing",
                    j = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            res <- c(list(x@sfilter), x@rfilters)
            names(res) <- sprintf("q=%i", seq(length(res) - 1, by = -1, length.out = length(res)))
            res <- res[j]
            if (length(res) == 1 & drop) {
              res <- res[[1]]
            }
            res
          })
#' @rdname filters_operations
#' @export
setMethod("[",
          signature(x = "finite_filters",
                    i = "ANY",
                    j = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            as.matrix(x)[i, j, ..., drop = drop]
          })
#' @export
to_seasonal.finite_filters <- function(x, s){
  x@sfilter <- to_seasonal(x@sfilter, s)
  x@rfilters <- unlist(lapply(x@rfilters, function(x){
    new_mm <- to_seasonal(x, s)
    rep(list(new_mm), s)
  }))
  x@lfilters <- unlist(lapply(x@lfilters, function(x){
    new_mm <- to_seasonal(x, s)
    rep(list(new_mm), s)
  }))
  x
}

#' Impute Incomplete Finite Filters
#'
#' @param x a [finite_filters()] object.
#' @param n integer specifying the number of imputed periods.
#' By default all the missing moving averages are imputed.
#' @param nperiod integer specifying how to imput missing date.
#' `nperiod = 1` means imputation using last filtered data (1 period backward),
#' `nperiod = 12` with monthly data means imputation using last year filtered data, etc.
#' @param backward,forward boolean indicating if the imputation should be done backward (on left filters), forward (on right filters).
#'
#' @details
#' When combining finite filters and a moving average, the first and/or the last points cannot be computed.
#'
#' For example, using the M2X12 moving average (symmetric moving average with coefficients \eqn{\theta = \begin{pmatrix} 1/24 & 1/12 & 1/12 & 1/12 & 1/12 & 1/12 & 1/12 & 1/12 & 1/12 & 1/12 & 1/12 & 1/12 & 1/24 \end{pmatrix}}), the first and last 6 points cannot be computed.
#'
#' `impute_last_obs()` allows to impute the first/last points using the `nperiod` previous filtered data. With `nperiod = 1`, the last filtered data is used for the imputation, with `nperiod = 12` and monthly data, the last year filtered data is used for the imputation, etc.
#'
#'
#' @examples
#' y <- window(retailsa$AllOtherGenMerchandiseStores, start = 2008)
#' M3 <- moving_average(rep(1/3, 3), lags = -1)
#' M3X3 <- M3 * M3
#' M2X12 <- (simple_ma(12, -6) + simple_ma(12, -5)) / 2
#' composite_ma <- M3X3 * M2X12
#' # The last 6 points cannot be computed
#' composite_ma
#' composite_ma * y
#' # they can be computed using the last filtered data
#' # e.g. to impute the first 3 missing months with last period:
#' impute_last_obs(composite_ma, n = 3, nperiod = 1) * y
#' # or using the filtered data of the same month in previous year
#' impute_last_obs(composite_ma, n = 6, nperiod = 12) * y
#' @export
impute_last_obs <- function(x, n, nperiod = 1, backward = TRUE, forward = TRUE) {
  if (is.moving_average(x))
    x <- finite_filters(sfilter = x)
  nrfilters <- length(x@rfilters)
  nlfilters <- length(x@lfilters)
  if (missing(n))
    n <- max(nrfilters, nlfilters)
  n_r <- min(upper_bound(x@sfilter) - nrfilters, n)
  n_l <- min(abs(lower_bound(x@sfilter)) - nlfilters, n)

  if (backward) {
    new_lfilters <- c(vector("list", n_l), x@lfilters)

    for (i in rev(seq_len(n_l))) {
      if (nperiod + i > length(new_lfilters)) {
        modified_filter <- x@sfilter
      } else {
        modified_filter <- new_lfilters[[i + nperiod]]
      }
      new_lfilters[[i]] <-
        modified_filter *
        moving_average(1, lags = nperiod)
    }
  } else {
    new_lfilters <- x@lfilters
  }

  if (forward) {
    new_rfilters <- c(x@rfilters, vector("list", n_r))
    for (i in seq_len(n_r)) {
      if (nrfilters + i - nperiod < 1) {
        modified_filter <- x@sfilter
      } else {
        modified_filter <- new_rfilters[[nrfilters + i - nperiod]]
      }
      new_rfilters[[nrfilters + i]] <-
        modified_filter *
        moving_average(1, lags = -nperiod)
    }
  } else {
    new_rfilters <- x@rfilters
  }
  finite_filters(x@sfilter, rfilters = new_rfilters, lfilters = new_lfilters)
}
