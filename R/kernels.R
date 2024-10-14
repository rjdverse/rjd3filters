#' Get the coefficients of a kernel
#'
#' Function to get the coefficient associated to a kernel. Those
#' coefficients are then used to compute the different filters.
#'
#' @inheritParams localpolynomials
#' @param sd_gauss standard deviation for gaussian kernel. By default 0.25.
#'
#' @return \code{tskernel} object (see \link[stats]{kernel}).
#' @export
#'
#' @examples
#' get_kernel("Henderson", horizon = 3)
get_kernel <- function(kernel = c("Henderson","Uniform", "Triangular",
                                  "Epanechnikov","Parabolic","BiWeight", "TriWeight","Tricube",
                                  "Trapezoidal", "Gaussian"),
                       horizon,
                       sd_gauss = 0.25){
  jkernel <- .r2jd_kernel(kernel, horizon, sd_gauss)
  coef <- sapply(as.integer(seq.int(from = 0, to = horizon, by = 1)),
                function(x) .jcall(jkernel, "D", "applyAsDouble", x))
  m <- horizon
  result <- list(coef = coef, m = m)
  attr(result, "name") <- kernel
  attr(result, "class") <- "tskernel"
  result
}
.r2jd_kernel <- function(
    kernel = c("Henderson","Uniform", "Triangular",
               "Epanechnikov","Parabolic","BiWeight", "TriWeight","Tricube",
               "Trapezoidal", "Gaussian"),
    horizon, sd_gauss = 0.25){

  if (is.null(kernel) || kernel[1]=="")
    return(.jnull("java/util/function/IntToDoubleFunction"))

  kernel <- match.arg(tolower(kernel)[1],
                      choices = c("henderson", "uniform", "triangular", "epanechnikov", "parabolic",
                                  "biweight", "triweight", "tricube", "trapezoidal", "gaussian"
                      ))
  if (kernel == "parabolic")
    kernel <- "epanechnikov"
  h <- as.integer(horizon)
  if (kernel == "gaussian"){
    jkernel <- .jcall("jdplus/toolkit/base/core/data/analysis/DiscreteKernel",
                      "Ljava/util/function/IntToDoubleFunction;",
                      tolower(kernel), h, sd_gauss)
  } else{
    jkernel <- .jcall("jdplus/toolkit/base/core/data/analysis/DiscreteKernel",
                      "Ljava/util/function/IntToDoubleFunction;",
                      tolower(kernel), h)
  }
  jkernel
}
