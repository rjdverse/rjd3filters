# Manipulation of moving averages

Manipulation of moving averages

## Usage

``` r
moving_average(
  x,
  lags = -length(x),
  trailing_zero = FALSE,
  leading_zero = FALSE
)

is.moving_average(x)

is_symmetric(x)

upper_bound(x)

lower_bound(x)

mirror(x)

# S3 method for class 'moving_average'
rev(x)

# S3 method for class 'moving_average'
length(x)

to_seasonal(x, s)

# S4 method for class 'moving_average'
show(object)
```

## Arguments

- x:

  vector of coefficients.

- lags:

  integer indicating the number of lags of the moving average.

- trailing_zero, leading_zero:

  boolean indicating whether to remove leading/trailing zero and NA.

- s:

  seasonal period for the `to_seasonal()` function.

- object:

  `moving_average` object.

## Details

A moving average is defined by a set of coefficient \\\boldsymbol
\theta=(\theta\_{-p},\dots,\theta\_{f})'\\ such all time series \\X_t\\
are transformed as: \$\$
M\_{\boldsymbol\theta}(X_t)=\sum\_{k=-p}^{+f}\theta_kX\_{t+k}=\left(\sum\_{k=-p}^{+f}\theta_kB^{-k}\right)X\_{t}
\$\$ The integer \\p\\ is defined by the parameter `lags`.

The function `to_seasonal()` transforms the moving average \\\boldsymbol
\theta\\ to: \$\$
M\_{\boldsymbol\theta'}(X_t)=\sum\_{k=-p}^{+f}\theta_kX\_{t+ks}=\left(\sum\_{k=-p}^{+f}\theta_kB^{-ks}\right)X\_{t}
\$\$

## Examples

``` r
y <- retailsa$AllOtherGenMerchandiseStores
e1 <- moving_average(rep(1,12), lags = -6)
e1 <- e1/sum(e1)
e2 <- moving_average(rep(1/12, 12), lags = -5)
M2X12 <- (e1 + e2)/2
stats::coef(M2X12)
#>        t-6        t-5        t-4        t-3        t-2        t-1          t 
#> 0.04166667 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 
#>        t+1        t+2        t+3        t+4        t+5        t+6 
#> 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.04166667 
M3 <- moving_average(rep(1/3, 3), lags = -1)
M3X3 <- M3 * M3
# M3X3 moving average applied to each month
M3X3
#> [1] "0.1111 B^2 + 0.2222 B + 0.3333 + 0.2222 F + 0.1111 F^2"
M3X3_seasonal <- to_seasonal(M3X3, 12)
# M3X3_seasonal moving average applied to the global series
M3X3_seasonal
#> [1] "0.1111 B^24 + 0.2222 B^12 + 0.3333 + 0.2222 F^12 + 0.1111 F^24"

def.par <- par(no.readonly = TRUE)
par(mai = c(0.5, 0.8, 0.3, 0))
layout(matrix(c(1,2), nrow = 1))
plot_gain(M3X3, main = "M3X3 applied to each month")
plot_gain(M3X3_seasonal, main = "M3X3 applied to the global series")

par(def.par)

# To apply the moving average
t <- y * M2X12
# Or use the filter() function:
t <- filter(y, M2X12)
si <- y - t
s <- si * M3X3_seasonal
# or equivalently:
s_mm <- M3X3_seasonal * (1 - M2X12)
s <- y * s_mm
graphics::plot(s)
```
