# Mean Square Revision Error (mmsre) filter

Provides an asymmetric filter based on the given reference filter
(usually symmetric) minimizing the mean square revision error.

## Usage

``` r
mmsre_filter(
  ref_filter,
  q,
  U,
  Z = NULL,
  delta = NULL,
  kernel = NULL,
  tweight = 0,
  passband = pi/12
)
```

## Arguments

- ref_filter:

  The reference filter (a
  [`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  object).

- q:

  The horizon of the asymmetric filter.

- U:

  Matrix of the constraints.

- Z:

  Matrix of the bias (can be `NULL`).

- delta:

  Coefficients of the linear model.

- kernel:

  The kernel used for weighting factors, by default, no weight is used.
  See
  [`lp_filter()`](https://rjdverse.github.io/rjd3filters/reference/lp_filter.md)
  for the available kernels.

- tweight:

  timeliness weight.

- passband:

  passband threshold.

## Details

The asymmetric filter \\\boldsymbol v=(v\_{-h},\dots,v{q})'\\ minimizes
the mean square revision error (mmsre) relative to the reference filter
\\\boldsymbol \theta=(\theta\_{-h},\dots,\theta\_{h'})'\\. The series
follows the model \$\$ \boldsymbol y=\boldsymbol U \boldsymbol \gamma +
\boldsymbol Z \boldsymbol \delta + \boldsymbol \varepsilon, \quad
\boldsymbol \varepsilon \sim \mathcal N(0,\sigma^2 \boldsymbol K^{-1}).
\$\$

With \\K\\ a set of weights (kernel), by default (`kernel = NULL`) no
weight is used. The matrix \\U\\ represents the constraints of the
symmetric filter (usually polynomials preservations), \\\boldsymbol
\theta\\, imposed to the asymmetric filter, \\\boldsymbol v\\.
Partitionning the matrix \\\boldsymbol U=\begin{pmatrix} \boldsymbol
U_p' & \boldsymbol U_f'\end{pmatrix}'\\ with \\\boldsymbol U_p\\ the
first \\h+q+1\\ rows and \\\boldsymbol U_f\\ the remaining rows, the
constraints are \\\boldsymbol U_p'\boldsymbol v=\boldsymbol
U'\boldsymbol \theta\\.

The matrix \\\boldsymbol Z\\ represents the bias of the asymmetric
filter: usually constraints imposed to the symmetric filter but not to
the asymmetric filter.

## References

Proietti, Tommaso and Alessandra Luati (2008). “Real time estimation in
local polynomial regression, with application to trend-cycle analysis”.

## See also

[`lp_filter()`](https://rjdverse.github.io/rjd3filters/reference/lp_filter.md).

## Examples

``` r
QL <- lp_filter(endpoints = "QL", ic = 3.5)
LC <- lp_filter(endpoints = "LC", ic = 3.5)
DAF <- lp_filter(endpoints = "DAF")
h6 <- QL[, "q=6"]
# To reproduce DAF filter
mmsre_filter(
  ref_filter = h6, q = 0,
  U = polynomial_matrix(l = - 6, d0 = 0, d1 = 3),
  kernel = "Henderson"
)
#> [1] " - 0.0172 B^6 + 0.0219 B^5 + 0.0400 B^4 - 0.0341 B^3 - 0.0979 B^2 + 0.1322 B + 0.9552"
DAF[, "q=0"]
#> [1] " - 0.0172 B^6 + 0.0219 B^5 + 0.0400 B^4 - 0.0341 B^3 - 0.0979 B^2 + 0.1322 B + 0.9552"
# To reproduce QL filter
mmsre_filter(
  ref_filter = h6, q = 1,
  delta = 2 / (sqrt(pi) * 3.5),
  U = polynomial_matrix(l = -6, d0 = 0, d1 = 1),
  Z = polynomial_matrix(l = -6, d0 = 2, d1 = 2)
)
#> [1] " - 0.0083 B^6 - 0.0395 B^5 - 0.0216 B^4 + 0.0466 B^3 + 0.1440 B^2 + 0.2392 B + 0.3058 + 0.3337 F"
QL[, "q=1"]
#> [1] " - 0.0083 B^6 - 0.0395 B^5 - 0.0216 B^4 + 0.0466 B^3 + 0.1440 B^2 + 0.2392 B + 0.3058 + 0.3337 F"

# Or using the Uniform kernel
mmsre_filter(
  ref_filter = h6, q = 2,
  # we multiply by the square root of the inverse of weights (1/13)
  # to get the same result as the QL filter
  delta = 2 / (sqrt(pi) * 3.5) * (sqrt(13)),
  U = polynomial_matrix(l = -6, d0 = 0, d1 = 0),
  Z = polynomial_matrix(l = -6, d0 = 1, d1 = 1),
  kernel = "Uniform"
)
#> [1] " - 0.0160 B^6 - 0.0249 B^5 + 0.0027 B^4 + 0.0678 B^3 + 0.1494 B^2 + 0.2160 B + 0.2414 + 0.2154 F + 0.1481 F^2"
LC[, "q=2"]
#> [1] " - 0.0160 B^6 - 0.0249 B^5 + 0.0027 B^4 + 0.0678 B^3 + 0.1494 B^2 + 0.2160 B + 0.2414 + 0.2154 F + 0.1481 F^2"
```
