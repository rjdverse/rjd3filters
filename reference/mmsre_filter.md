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
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): RcallMethod: cannot determine object class
LC <- lp_filter(endpoints = "LC", ic = 3.5)
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
DAF <- lp_filter(endpoints = "DAF")
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): RcallMethod: cannot determine object class
h6 <- QL[, "q=6"]
#> Error: object 'QL' not found
# To reproduce DAF filter
mmsre_filter(
  ref_filter = h6, q = 0,
  U = polynomial_matrix(l = - 6, d0 = 0, d1 = 3),
  kernel = "Henderson"
)
#> Error: object 'h6' not found
DAF[, "q=0"]
#> Error: object 'DAF' not found
# To reproduce QL filter
mmsre_filter(
  ref_filter = h6, q = 1,
  delta = 2 / (sqrt(pi) * 3.5),
  U = polynomial_matrix(l = -6, d0 = 0, d1 = 1),
  Z = polynomial_matrix(l = -6, d0 = 2, d1 = 2)
)
#> Error: object 'h6' not found
QL[, "q=1"]
#> Error: object 'QL' not found

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
#> Error: object 'h6' not found
LC[, "q=2"]
#> Error: object 'LC' not found
```
