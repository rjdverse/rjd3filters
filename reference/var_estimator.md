# Variance Estimator

Variance Estimator

## Usage

``` r
var_estimator(x, coef, ...)
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

## Details

Let \\(\theta_i)\_{-p\leq i \leq q}\\ be a moving average of length
\\p+q+1\\ used to filter a time series \\(y_i)\_{1\leq i \leq n}\\. It
is equivalent to a local regression and the associated error variance
\\\sigma^2\\ can be estimated using the normalized residual sum of
squares, which can be simplified as: \$\$
\hat\sigma^2=\frac{1}{n-(p+q)}\sum\_{t=p+1}^{n-q} \frac{(y_t-\hat
\mu_t)^2}{1-2w_0^2+\sum\_{i=-p}^q w_i^2} \$\$

## References

Loader, Clive. 1999. Local regression and likelihood. New York:
Springer-Verlag.

## See also

[`df_var()`](https://rjdverse.github.io/rjd3filters/reference/df_var.md).
