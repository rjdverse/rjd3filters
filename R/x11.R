#' X-11 Decomposition With Custom Trend Filters
#'
#' Perform the X-11 decomposition using custom trend filter
#' @param x input time-series.
#' @param period period.
#' @param trend.coefs coefficients of the filters used for the trend-cycle extraction from
#' the real-time asymmetric filter to the symmetric filter. Can be a, object of class \code{"list"},
#' \code{"matrix"}, \code{"lp_filter"} or \code{"rkhs_filter"}.
#' @param mul boolean indicating if the decomposition mode is multiplicative.
#' @param seas.s0,seas.s1 seasonal filters.
#' @param extreme.lsig,extreme.usig boundaries used for outlier correction in irregular.
#' @param userdefined a vector containing the additional output variables.
#' @examples
#' x <- retailsa$AllOtherGenMerchandiseStores
#' decomposition_lp <- x11(x, trend.coefs = lp_filter())
#' decomposition_rkhs <- x11(x, trend.coefs = rkhs_filter())
#' plot(x)
#' lines(decomposition_lp$decomposition[,"t"], col = "red")
#' lines(decomposition_rkhs$decomposition[,"t"], col = "green")
#' @importFrom stats ts.union
#' @export
x11 <- function(x, period = frequency(x),
                trend.coefs,  mul=TRUE,
                seas.s0=c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"),
                seas.s1=c("S3X5", "S3X3", "S3X1", "S3X9", "S3X15"),
                extreme.lsig=1.5, extreme.usig=2.5, userdefined = NULL){
  seas.s0<-match.arg(toupper(seas.s0)[1], choices = c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"))
  seas.s1<-match.arg(toupper(seas.s1)[1], choices = c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"))

  if (!inherits(trend.coefs, "finite_filters")) {
    trend.coefs <- finite_filters(trend.coefs)
  }

  sym_filter <- trend.coefs@sfilter
  asy_filter <- trend.coefs@lfilters
  leftTrendFilter <- do.call(cbind, lapply(asy_filter, function(x){
    c(coef(x), rep(0, 2*length(asy_filter)- length(x)))
  }))
  leftTrendFilter <- leftTrendFilter[,seq(ncol(leftTrendFilter), 1)]
  if(length(sym_filter) != 2*ncol(leftTrendFilter)+1){
    stop(sprintf("The symmetric filter is of length %i but only %i asymmetric filters provided",
                 length(sym_filter),
                 2*ncol(leftTrendFilter)+1))
  }
  ltrendf <- rjd3toolkit::.r2jd_matrix(leftTrendFilter)

  ctrendf <- .jcall("jdplus/toolkit/base/api/data/DoubleSeq",
                   "Ljdplus/toolkit/base/api/data/DoubleSeq;",
                   "of", coef(sym_filter))
  # # check :
  # jf <- J("jdplus.filters.base.core.filters.Filtering")$of(ctrendf, ltrendf)
  # .jd2r_finitefilters(jf)
  # trend.coefs

  x11decomp <- J("jdplus/x11plus/base/r/X11Decomposition")

  jrslt <- x11decomp$trendX11(as.numeric(x), period, mul,
                             ctrendf, ltrendf,
                             seas.s0, seas.s1, extreme.lsig, extreme.usig)
  decomposition <- cbind(x,
                         rjd3toolkit::.proc_vector(jrslt, "d11"),
                         rjd3toolkit::.proc_vector(jrslt, "d12"),
                         rjd3toolkit::.proc_vector(jrslt, "d10"),
                         rjd3toolkit::.proc_vector(jrslt, "d13"))
  colnames(decomposition) <- c("y", "sa", "t", "s", "i")
  res <- structure(list(
    decomposition=decomposition),
    class="X11_RSLT")

  res$user_defined <- rjd3toolkit::user_defined(rjd3toolkit::.jd3_object(jrslt, result = TRUE),
                                               userdefined = userdefined)
  return(res)
}
