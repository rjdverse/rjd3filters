# FST criteria

Compute the Fidelity, Smoothness and Timeliness (FST) criteria

## Usage

``` r
fst(weights, lags, passband = pi/6, ...)
```

## Arguments

- weights:

  either a `"moving_average"` or a numeric vector containing weights.

- lags:

  Lags of the moving average (when `weights` is not a
  `"moving_average"`).

- passband:

  Passband threshold for timeliness criterion.

- ...:

  other unused arguments.

## Value

The values of the 3 criteria, the gain and phase of the associated
filter.

## References

Grun-Rehomme, Michel, Fabien Guggemos, and Dominique Ladiray (2018).
“Asymmetric Moving Averages Minimizing Phase Shift”. In: Handbook on
Seasonal Adjustment,
<https://ec.europa.eu/eurostat/web/products-manuals-and-guidelines/-/ks-gq-18-001>.

## Examples

``` r
filter <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "LC")
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): RcallMethod: cannot determine object class
fst(filter[, "q=0"])
#> Error in filter[, "q=0"]: object of type 'closure' is not subsettable
# To compute the statistics on all filters:
fst(filter)
#> Error in .jcall("jdplus/filters/base/core/AdvancedFiltersToolkit", "Ljdplus/filters/base/core/AdvancedFiltersToolkit$FSTResult;",     "fst", weights, as.integer(lags), passband): java.lang.UnsupportedClassVersionError: jdplus/filters/base/r/LocalPolynomialFilters has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
```
