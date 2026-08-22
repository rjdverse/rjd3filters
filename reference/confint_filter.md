# Confidence intervals

Confidence intervals

## Usage

``` r
confint_filter(
  x,
  coef,
  coef_var = coef,
  level = 0.95,
  asymmetric_var = TRUE,
  gaussian_distribution = FALSE,
  exact_df = TRUE,
  ...
)
```

## Arguments

- x:

  input time series.

- coef:

  moving-average
  ([`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md))
  or finite filter
  ([`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md))
  used to filter the series.

- coef_var:

  moving-average
  ([`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md))
  or finite filter
  ([`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md))
  used compute the variance (throw
  [`var_estimator()`](https://rjdverse.github.io/rjd3filters/reference/var_estimator.md)).
  By default equal to `coef`.

- level:

  confidence level.

- asymmetric_var:

  if `asymmetric_var = TRUE` then the variance is estimated for each
  asymmetric filters instead of using the variance associated the
  symmetric estimates.

- gaussian_distribution:

  if `TRUE` use the normal distribution to compute the confidence
  interval, otherwise use the t-distribution.

- exact_df:

  if `TRUE` compute the exact degrees of freedom for the t-distribution
  (when `gaussian_distribution = FALSE`), otherwise uses an
  approximation.

- ...:

  other arguments passed to the function
  [`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  to convert `coef` to a `"moving_average"` object.

## Details

Let \\(\theta_i)\_{-p\leq i \leq q}\\ be a moving average of length
\\p+q+1\\ used to filter a time series \\(y_i)\_{1\leq i \leq n}\\. Let
denote \\\hat{\mu}\_t\\ the filtered series computed at time \\t\\ as:
\$\$ \hat{\mu}\_t = \sum\_{i=-p}^q \theta_i y\_{t+i}. \$\$ If
\\\hat{\mu}\_t\\ is unbiased, a approximate confidence for the true mean
is: \$\$ \left\[\hat{\mu}\_t - z\_{1-\alpha/2} \hat{\sigma}
\sqrt{\sum\_{i=-p}^q\theta_i^2}; \hat{\mu}\_t + z\_{1-\alpha/2}
\hat{\sigma} \sqrt{\sum\_{i=-p}^q\theta_i^2} \right\], \$\$ where
\\z\_{1-\alpha/2}\\ is the quantile \\1-\alpha/2\\ of the standard
normal distribution.

The estimate of the variance \\\hat{\sigma}\\ is obtained using
[`var_estimator()`](https://rjdverse.github.io/rjd3filters/reference/var_estimator.md)
with the parameter `coef_var`. The assumption that \\\hat{\mu}\_t\\ is
unbiased is rarely exactly true, so variance estimates and confidence
intervals are usually computed at small bandwidths where bias is small.

When `coef` (or `coef_var`) is a finite filter, the last points of the
confidence interval are computed using the corresponding asymmetric
filters

## References

Loader, Clive. 1999. Local regression and likelihood. New York:
Springer-Verlag.

## Examples

``` r
x <- retailsa$DrinkingPlaces
coef <- lp_filter(6)
confint <- confint_filter(x, coef)
graphics::plot(confint, plot.type = "single",
     col = c("red", "black", "black"),
     lty = c(1, 2, 2))
```
