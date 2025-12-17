# Diagnostics and goodness of fit of filtered series

Set of functions to compute diagnostics and goodness of fit of filtered
series: cross validation (`cv()`) and cross validate estimate (`cve()`),
leave-one-out cross validation estimate (`loocve`), CP statistic
(`cp()`) and Rice's T statistics (`rt()`).

## Usage

``` r
cve(x, coef, ...)

cv(x, coef, ...)

loocve(x, coef, ...)

rt(x, coef, ...)

cp(x, coef, var, ...)
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

- var:

  variance used to compute the CP statistic (`cp()`).

## Details

Let \\(\theta_i)\_{-p\leq i \leq q}\\ be a moving average of length
\\p+q+1\\ used to filter a time series \\(y_i)\_{1\leq i \leq n}\\. Let
denote \\\hat{\mu}\_t\\ the filtered series computed at time \\t\\ as:
\$\$ \hat{\mu}\_t = \sum\_{i=-p}^q \theta_i y\_{t+i}. \$\$

The cross validation estimate (`cve()`) is defined as the time series
\\Y_t-\hat{\mu}\_{-t}\\ where \\\hat{\mu}\_{-t}\\ is the leave-one-out
cross validation estimate (`loocve()`) defined as the filtered series
computed deleting the observation \\t\\ and remaining all the other
points. The cross validation statistics (`cv()`) is defined as: \$\$
CV=\frac{1}{n-(p+q)} \sum\_{t=p+1}^{n-q} \left(y_t -
\hat{\mu}\_{-t}\right)^2. \$\$ In the case of filtering with a moving
average, we can show that: \$\$ \hat{\mu}\_{-t}= \frac{\hat{\mu}\_t -
\theta_0 y_t}{1-\theta_0} \$\$ and \$\$ CV=\frac{1}{n-(p+q)}
\sum\_{t=p+1}^{n-q} \left(\frac{y_t -
\hat{\mu}\_{t}}{1-\theta_0}\right)^2. \$\$

In the case of filtering with a moving average, the CP estimate of risk
(introduced by Mallows (1973); `cp()`) can be defined as: \$\$
CP=\frac{1}{\sigma^2} \sum\_{t=p+1}^{n-q} \left(y_t -
\hat{\mu}\_{t}\right)^2 -(n-(p+q))(1-2\theta_0). \$\$ The CP method
requires an estimate of \\\sigma^2\\ (`var` parameter). The usual use of
CP is to compare several different fits (for example different
bandwidths): one should use the same estimate of \\\hat{\sigma}^2\\ for
all fits (using for example
[`var_estimator()`](https://rjdverse.github.io/rjd3filters/reference/var_estimator.md)).
The recommendation of Cleveland and Devlin (1988) is to compute
\\\hat{\sigma}^2\\ from a fit at the smallest bandwidth under
consideration, at which one should be willing to assume that bias is
negligible.

The Rice's T statistic (`rt()`) is defined as: \$\$ \frac{1}{n-(p+q)}
\sum\_{t=p+1}^{n-q} \frac{ \left(y_t - \hat{\mu}\_{t}\right)^2 }{
1-2\theta_0 } \$\$

## References

Loader, Clive. 1999. Local regression and likelihood. New York:
Springer-Verlag.

Mallows, C. L. (1973). Some comments on Cp. Technometrics 15, 661– 675.

Cleveland, W. S. and S. J. Devlin (1988). Locally weighted regression:
An approach to regression analysis by local fitting. Journal of the
American Statistical Association 83, 596–610.
