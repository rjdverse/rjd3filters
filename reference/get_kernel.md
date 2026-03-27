# Get the coefficients of a kernel

Function to get the coefficient associated to a kernel. Those
coefficients are then used to compute the different filters.

## Usage

``` r
get_kernel(
  kernel = c("Henderson", "Uniform", "Triangular", "Epanechnikov", "Parabolic",
    "BiWeight", "TriWeight", "Tricube", "Trapezoidal", "Gaussian"),
  horizon,
  sd_gauss = 0.25
)
```

## Arguments

- kernel:

  kernel uses.

- horizon:

  horizon (bandwidth) of the symmetric filter.

- sd_gauss:

  standard deviation for gaussian kernel. By default 0.25.

## Value

`tskernel` object (see [kernel](https://rdrr.io/r/stats/kernel.html)).

## Examples

``` r
get_kernel("Henderson", horizon = 3)
#> Error in .jcall("jdplus/toolkit/base/core/data/analysis/DiscreteKernel",     "Ljava/util/function/IntToDoubleFunction;", tolower(kernel),     h): java.lang.UnsupportedClassVersionError: jdplus/filters/base/core/AdvancedFiltersToolkit has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
```
