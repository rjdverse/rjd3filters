# Apply Local Polynomials Filters

Apply Local Polynomials Filters

## Usage

``` r
localpolynomials(
  x,
  horizon = 6,
  degree = 3,
  kernel = c("Henderson", "Uniform", "Biweight", "Trapezoidal", "Triweight", "Tricube",
    "Gaussian", "Triangular", "Parabolic"),
  endpoints = c("LC", "QL", "CQ", "CC", "DAF"),
  ic = 3.5,
  tweight = 0,
  passband = pi/12
)
```

## Arguments

- x:

  input time-series.

- horizon:

  horizon (bandwidth) of the symmetric filter.

- degree:

  degree of polynomial.

- kernel:

  kernel uses.

- endpoints:

  method for endpoints.

- ic:

  ic ratio.

- tweight:

  timeliness weight.

- passband:

  passband threshold.

## Value

the target signal

## References

Proietti, Tommaso and Alessandra Luati (2008). “Real time estimation in
local polynomial regression, with application to trend-cycle analysis”.

## See also

[`lp_filter()`](https://rjdverse.github.io/rjd3filters/reference/lp_filter.md).

## Examples

``` r
x <- retailsa$AllOtherGenMerchandiseStores
trend <- localpolynomials(x, horizon = 6)
graphics::plot(x)
graphics::lines(trend, col = "red")
```
