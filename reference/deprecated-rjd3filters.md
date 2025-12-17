# Deprecated function

Deprecated function

## Usage

``` r
cross_validation(x, coef, ...)

implicit_forecast(x, coefs)
```

## Arguments

- x:

  input time series.

- coef:

  vector of coefficients or a moving-average
  ([`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)).

- ...:

  other arguments passed to the function
  [`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  to convert `coef` to a `"moving_average"` object.
