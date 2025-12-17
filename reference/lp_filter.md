# Local Polynomials Filters

Local Polynomials Filters

## Usage

``` r
lp_filter(
  horizon = 6,
  degree = 3,
  kernel = c("Henderson", "Uniform", "Biweight", "Trapezoidal", "Triweight", "Tricube",
    "Gaussian", "Triangular", "Parabolic"),
  endpoints = c("LC", "QL", "CQ", "CC", "DAF", "CN"),
  ic = 3.5,
  tweight = 0,
  passband = pi/12
)
```

## Arguments

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

a
[`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md)
object.

## Details

- "LC": Linear-Constant filter

- "QL": Quadratic-Linear filter

- "CQ": Cubic-Quadratic filter

- "CC": Constant-Constant filter

- "DAF": Direct Asymmetric filter

- "CN": Cut and Normalized Filter

## References

Proietti, Tommaso and Alessandra Luati (2008). “Real time estimation in
local polynomial regression, with application to trend-cycle analysis”.

## See also

[`mmsre_filter()`](https://rjdverse.github.io/rjd3filters/reference/mmsre_filter.md)
[`localpolynomials()`](https://rjdverse.github.io/rjd3filters/reference/localpolynomials.md).

## Examples

``` r
henderson_f <- lp_filter(horizon = 6, kernel = "Henderson")
plot_coef(henderson_f)
```
