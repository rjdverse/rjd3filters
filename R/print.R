#' @rdname moving_average
#' @importFrom methods show
#' @export
setMethod(
    f = "show",
    signature = "moving_average",
    definition = function(object) {
        print(.jcall(.ma2jd(object), "S", "toString"))
        invisible(object)
    }
)


#' Print method for Finite Filters Objects
#'
#' Printing objects of class `"finite_filters"`.
#'
#' @param x,object An object of class `"finite_filters"`.
#' @param sfilter,rfilters,lfilters Logical indicating whether to add the symmetric filter
#'   (`sfilter = TRUE` by default), the right filters (`rfilters = TRUE` by default) or
#'   the left filters (`lfilters = FALSE` by default).
#' @param zero_as_na Logical. Asymmetric filters have fewer terms than the
#'   symmetric filter. If `TRUE`, the missing past/future terms in the
#'   matrix representation are filled with `NA` instead of `0`.
#'   Defaults to `FALSE`.
#' @param ... Further arguments to be passed to or from methods.
#'
#' @details
#' The print method just calls the `as.matrix()` method and prints the resulting matrix.
#'
#' @returns The print method returns the argument `x`, invisibly, as for all [base::print()] methods.
#'
#' `as.matrix.finite_filters()` returns a matrix containing the coefficients of the symmetric and asymmetric filters.
#' The columns are named according to the number of future points used (for right filters).
#' For example, column `q=0` corresponds to the asymmetric filter using no future points,
#' while `q=1` corresponds to the filter using one future point, and so on.
#' A similar convention with negative values is used for left filters (when printed).
#'
#' @aliases print,finite_filters-method
#' @rdname print.finite_filters
#' @seealso [finite_filters()].
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' ff_lp <- lp_filter()
#' as.matrix(ff_lp)
#' @export
setMethod(
    f = "print",
    signature = "finite_filters",
    definition = function(x, sfilter = TRUE,
                          rfilters = TRUE,
                          lfilters = FALSE,
                          zero_as_na = FALSE, ...) {
        m <- as.matrix(x, sfilter = sfilter,
                       rfilters = rfilters,
                       lfilters = lfilters,
                       zero_as_na = zero_as_na)
        print(m, ...)
        invisible(x)
    }
)
#' @importFrom methods show
#' @aliases show,finite_filters-method
#' @rdname print.finite_filters
#' @export
setMethod(
    f = "show",
    signature = "finite_filters",
    definition = function(object) {
        x <- as.matrix(object)
        print(x)
    }
)

#' @method as.matrix finite_filters
#' @rdname print.finite_filters
#' @export
as.matrix.finite_filters <- function(
        x,
        sfilter = TRUE,
        rfilters = TRUE,
        lfilters = FALSE,
        zero_as_na = FALSE,
        ...
) {
    sfilter_s <- rfilters_s <- lfilters_s <-
        index_s <- index_r <- index_l <- NULL
    if (!any(sfilter, rfilters, lfilters)) {
        return(NULL)
    }
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
    mat <- do.call(
        cbind,
        c(lfilters_s, sfilter_s, rfilters_s, list(zero_as_na = zero_as_na))
    )
    colnames(mat) <- sprintf("q=%i", c(index_l, index_s, index_r))
    mat
}
