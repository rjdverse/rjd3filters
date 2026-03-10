# Accuracy/smoothness/timeliness criteria through spectral decomposition

Accuracy/smoothness/timeliness criteria through spectral decomposition

## Usage

``` r
mse(aweights, sweights, density = c("uniform", "rw"), passband = pi/6, ...)
```

## Arguments

- aweights:

  `moving_average` object or weights of the asymmetric filter (from -n
  to m).

- sweights:

  `moving_average` object or weights of the symmetric filter (from 0 to
  n or -n to n).

- density:

  hypothesis on the spectral density: `"uniform"` (= white woise, the
  default) or `"rw"` (= random walk).

- passband:

  passband threshold.

- ...:

  other unused arguments.

## Value

The criteria

## References

Wildi, Marc and McElroy, Tucker (2019). “The trilemma between accuracy,
timeliness and smoothness in real-time signal extraction”. In:
International Journal of Forecasting 35.3, pp. 1072–1084.

## Examples

``` r
filter <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "LC")
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): RcallMethod: cannot determine object class
sweights <- filter[, "q=6"]
#> Error in filter[, "q=6"]: object of type 'closure' is not subsettable
aweights <- filter[, "q=0"]
#> Error in filter[, "q=0"]: object of type 'closure' is not subsettable
mse(aweights, sweights)
#> Error: object 'aweights' not found
# Or to compute directly the criteria on all asymmetric filters:
mse(filter)
#> Error in mse.default(filter): argument "sweights" is missing, with no default
```
