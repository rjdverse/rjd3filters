# Retrieve underlying forecasts corresponding to the asymmetric filters

Function to retrieve the underlying forecasts corresponding to the
asymmetric filters: the asymmetric filtering produce the same result
than the symmetric filters applied to the extended series.

## Usage

``` r
underlying_forecasts(x, coefs)
```

## Arguments

- x:

  a univariate or multivariate time series.

- coefs:

  a `matrix` or a `list` that contains all the coefficients of the
  asymmetric and symmetric filters. (from the symmetric filter to the
  shortest). See details.

## Details

Let \\h\\ be the bandwidth of the symmetric filter, \\v\_{-h}, \ldots,
v_h\\ the coefficients of the symmetric filter and \\w\_{-h}^q, \ldots,
w_h^q\\ the coefficients of the asymmetric filter used to estimate the
trend when \\q\\ future values are known (with the convention
\\w\_{q+1}^q=\ldots=w_h^q=0\\). Let denote \\y\_{-2h+1},\ldots, y_0\\
the last \\2h\\ available values of the input times series. The
underlying forecast \\y\_{1}^\*,\dots y_h^\*\\ induced by \\w^0,\dots
w^{h-1}\\ are defined by: \$\$ \forall q\in\\0,...,h-1\\, \quad
\sum\_{i=-h}^q v_iy\_{i-q} + \sum\_{i=q+1}^h v_iy\_{i-q}^\*
=\sum\_{i=-h}^q w_i^qy\_{i-q}. \$\$ Note that this is solved
numerically: the solution isn't exact.

## Examples

``` r
x <- retailsa$AllOtherGenMerchandiseStores

ql <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "QL")
lc <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "LC")
f_ql <- underlying_forecasts(x, ql)
f_lc <- underlying_forecasts(x, lc)

graphics::plot(window(x, start = 2007), xlim = c(2007, 2012))
graphics::lines(
    stats::ts(
        c(utils::tail(x, 1), f_ql),
        frequency = stats::frequency(x),
        start = stats::end(x)
    ),
    col = "red",
    lty = 2
)
graphics::lines(
    stats::ts(
        c(utils::tail(x, 1), f_lc),
        frequency = stats::frequency(x),
        start = stats::end(x)
    ),
    col = "blue",
    lty = 2
)
```
