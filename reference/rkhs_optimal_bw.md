# Optimal Bandwith of Reproducing Kernel Hilbert Space (RKHS) Filters

Function to export the optimal bandwidths used in Reproducing Kernel
Hilbert Space (RKHS) filters

## Usage

``` r
rkhs_optimal_bw(
  horizon = 6,
  degree = 2,
  kernel = c("Biweight", "Henderson", "Epanechnikov", "Triangular", "Uniform",
    "Triweight"),
  asymmetricCriterion = c("Timeliness", "FrequencyResponse", "Accuracy", "Smoothness"),
  density = c("uniform", "rw"),
  passband = 2 * pi/12,
  optimal.minBandwidth = horizon,
  optimal.maxBandwidth = 3 * horizon
)
```

## Arguments

- horizon:

  horizon (bandwidth) of the symmetric filter.

- degree:

  degree of polynomial.

- kernel:

  kernel uses.

- asymmetricCriterion:

  the criteria used to compute the optimal bandwidth. If `"Undefined"`,
  \\m+1\\ is used.

- density:

  hypothesis on the spectral density: `"uniform"` (= white woise, the
  default) or `"rw"` (= random walk).

- passband:

  passband threshold.

- optimal.minBandwidth, optimal.maxBandwidth:

  the range used for the optimal bandwith selection.

## Examples

``` r
rkhs_optimal_bw(asymmetricCriterion = "Timeliness")
#> Error in .jcall("jdplus/filters/base/r/RKHSFilters", "[D", "optimalBandwidth",     as.integer(horizon), as.integer(degree), kernel, asymmetricCriterion,     density == "rw", passband, optimal.minBandwidth, optimal.maxBandwidth): RcallMethod: cannot determine object class
rkhs_optimal_bw(asymmetricCriterion = "Timeliness", optimal.minBandwidth = 6.2)
#> Error in .jcall("jdplus/filters/base/r/RKHSFilters", "[D", "optimalBandwidth",     as.integer(horizon), as.integer(degree), kernel, asymmetricCriterion,     density == "rw", passband, optimal.minBandwidth, optimal.maxBandwidth): java.lang.UnsupportedClassVersionError: jdplus/filters/base/r/RKHSFilters has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
```
