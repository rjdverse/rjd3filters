# Compute the degrees of freedom for the variance estimator

Compute the degrees of freedom for the variance estimator

## Usage

``` r
df_var(n, coef, exact_df = FALSE)
```

## Arguments

- n:

  number of observations

- coef:

  moving average
  ([`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md))
  used to filter the series.

- exact_df:

  if `TRUE` compute the exact degrees of freedom for the t-distribution
  (when `gaussian_distribution = FALSE`), otherwise uses an
  approximation.

## See also

[`var_estimator()`](https://rjdverse.github.io/rjd3filters/reference/var_estimator.md).
