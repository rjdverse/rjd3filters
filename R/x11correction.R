#' 'X11' Extreme Values Corrector
#'
#' @param x the analysed time series.
#' @param period the period of the input time series if `x` is not a `"ts"` object.
#' @param corrected_s other time series if the series being corrected is different from x.
#' @param lsigma,usigma the lower and the upper sigma boundaries for the detection of extreme values.
#' @param mul boolean indicating if the decomposition is multiplicative or additive.
#' @param start position of the first "complete" considered period.
#' @param clean_extremities boolean indicating if the extremities should be cleaned.
#'
#' @details
#' The 'X11' Extreme Values Corrector is used to compute the tables
#' b4, b4g, b9, b9g, b17, b20, c17 and c20.
#'
#' @importFrom stats end window
#' @export
x11_extreme_values_corr <- function(x,
                                    corrected_s,
                                    period,
                                    lsigma = 1.5,
                                    usigma = 2.5,
                                    mul = FALSE,
                                    start = 0,
                                    clean_extremities = TRUE) {
  if (missing(period)) {
    if (is.ts(x)) {
      period <- stats::frequency(x)
    } else {
      stop("The period parameter must be defined")
    }
  }

  P <- .jcast(new( J("java.lang.Double"), as.character(period)),
              "java.lang.Number")
  dmode <- .jcall("jdplus/sa/base/api/DecompositionMode",
                  "Ljdplus/sa/base/api/DecompositionMode;",
                  "valueOf",
                  ifelse(mul, "mul", "Additive"))
  if (clean_extremities) {
    x_c <- rjd3toolkit::clean_extremities(x)
  } else {
    x_c  <- x
  }

  jx <- .r2jd_doubleseq(x_c)
  x11context = J("jdplus/x11plus/base/core/X11Context")$
    builder()$
    mode(dmode)$
    period(P)$
    lowerSigma(lsigma)$
    upperSigma(usigma)$
    build()
  dfvc <- .jnew("jdplus/x11plus/base/core/DefaultExtremeValuesCorrector")
  .jcall(dfvc, "V",
         "setStart", as.integer(start))
  .jcall(dfvc, "V",
         "analyse", jx, x11context)
  .jcall(dfvc, "V",
         "analyse", jx, x11context)
  obs_w <- .jcall(.jcall(dfvc, "Ljdplus/toolkit/base/api/data/DoubleSeq;",
                         "getObservationWeights"), "[D", "toArray")
  if (missing(corrected_s)) {
    corr_f <-  .jcall(.jcall(dfvc, "Ljdplus/toolkit/base/api/data/DoubleSeq;",
                             "getCorrectionFactors"), "[D", "toArray")
  } else {
    if (clean_extremities) {
      corrected_s <- rjd3toolkit::clean_extremities(corrected_s)
    }
    corr_f <- .jcall(.jcall(dfvc, "Ljdplus/toolkit/base/api/data/DoubleSeq;",
                            "computeCorrections",
                            .r2jd_doubleseq(corrected_s)), "[D", "toArray")
  }
  res <- cbind(obs_w,
               corr_f
  )
  colnames(res) <- c("obs_weight", "correction_factors")
  if (is.ts(x)) {
    res <- ts(res, start = start(x_c), frequency = frequency(x_c))
    res <- window(res, start = start(x), end = end(x), extend = TRUE)
  }
  res
}
