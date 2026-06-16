# Manipulating Finite Filters

Manipulating Finite Filters

## Usage

``` r
finite_filters(
  sfilter,
  rfilters = NULL,
  lfilters = NULL,
  first_to_last = FALSE
)

is.finite_filters(x)

# S4 method for class 'finite_filters'
show(object)
```

## Arguments

- sfilter:

  the symmetric filter
  ([`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  object) or a `matrix` or `list` with all the coefficients.

- rfilters:

  the right filters (used on the last points).

- lfilters:

  the left filters (used on the first points).

- first_to_last:

  boolean indicating if the first element of `rfilters` is the first
  asymmetric filter (when only one observation is missing) or the last
  one (real-time estimates).

- x:

  object to test the class.

- object:

  `finite_filters` object.

## Examples

``` r
ff_lp <- lp_filter()
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): RcallMethod: cannot determine object class
ff_simple_ma <- finite_filters(moving_average(c(1, 1, 1), lags = -1)/3,
               rfilters = list(moving_average(c(1, 1), lags = -1)/2))
#> Error in .jfindClass(as.character(class), class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
ff_lp
#> Error: object 'ff_lp' not found
ff_simple_ma
#> Error: object 'ff_simple_ma' not found
ff_lp * ff_simple_ma
#> Error: object 'ff_lp' not found
```
