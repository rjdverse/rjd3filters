# Direct Filter Approach

Direct Filter Approach

## Usage

``` r
dfa_filter(
  horizon = 6,
  degree = 0,
  density = c("uniform", "rw"),
  targetfilter = lp_filter(horizon = horizon)@sfilter,
  passband = 2 * pi/12,
  accuracy.weight = 1/3,
  smoothness.weight = 1/3,
  timeliness.weight = 1/3
)
```

## Arguments

- horizon:

  horizon (bandwidth) of the symmetric filter.

- degree:

  degree of polynomial.

- density:

  hypothesis on the spectral density: `"uniform"` (= white woise, the
  default) or `"rw"` (= random walk).

- targetfilter:

  the weights of the symmetric target filters (by default the Henderson
  filter).

- passband:

  passband threshold.

- accuracy.weight, smoothness.weight, timeliness.weight:

  the weight used for the optimisation. The weight associated to the
  residual is derived so that the sum of the four weights are equal to
  1.

## Details

Moving average computed by a minimisation of a weighted mean of three
criteria under polynomials constraints. The criteria come from the
decomposition of the mean squared error between th trend-cycle

Let \\\boldsymbol \theta=(\theta\_{-p},\dots,\theta\_{f})'\\ be a moving
average where \\p\\ and \\f\\ are two integers defined by the parameter
`lags` and `leads`. The three criteria are:

## Examples

``` r
dfa_filter(horizon = 6, degree = 0)
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
dfa_filter(horizon = 6, degree = 2)
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): RcallMethod: cannot determine object class
```
