# Optimal Bandwidth of Reproducing Kernel Hilbert Space (RKHS) Filters

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

  hypothesis on the spectral density: `"uniform"` (= white noise, the
  default) or `"rw"` (= random walk).

- passband:

  passband threshold.

- optimal.minBandwidth, optimal.maxBandwidth:

  the range used for the optimal bandwidth selection.

## Examples

``` r
rkhs_optimal_bw(asymmetricCriterion = "Timeliness")
#>    q=0    q=1    q=2    q=3    q=4    q=5 
#> 6.0000 6.0000 6.3875 8.1500 9.3500 6.0000 
rkhs_optimal_bw(asymmetricCriterion = "Timeliness", optimal.minBandwidth = 6.2)
#>       q=0       q=1       q=2       q=3       q=4       q=5 
#>  6.200000  6.200000  6.384375  8.148229  9.352812 10.391458 
```
