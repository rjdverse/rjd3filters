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

  boolean indicating wheter to remove leading/trailing zero and NA.

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
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
e2 <- moving_average(rep(1/12, 12), lags = -5)
M2X12 <- (e1 + e2)/2
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
coef(M2X12)
#> Error: object 'M2X12' not found
M3 <- moving_average(rep(1/3, 3), lags = -1)
M3X3 <- M3 * M3
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
# M3X3 moving average applied to each month
M3X3
#> Error: object 'M3X3' not found
M3X3_seasonal <- to_seasonal(M3X3, 12)
#> Error: object 'M3X3' not found
# M3X3_seasonal moving average applied to the global series
M3X3_seasonal
#> Error: object 'M3X3_seasonal' not found

def.par <- par(no.readonly = TRUE)
par(mai = c(0.5, 0.8, 0.3, 0))
layout(matrix(c(1,2), nrow = 1))
plot_gain(M3X3, main = "M3X3 applied to each month")
#> Error: object 'M3X3' not found
plot_gain(M3X3_seasonal, main = "M3X3 applied to the global series")
#> Error: object 'M3X3_seasonal' not found
par(def.par)

# To apply the moving average
t <- y * M2X12
#> Error: object 'M2X12' not found
# Or use the filter() function:
t <- filter(y, M2X12)
#> Error: object 'M2X12' not found
si <- y - t
#> Error in `-.default`(y, t): non-numeric argument to binary operator
s <- si * M3X3_seasonal
#> Error: object 'si' not found
# or equivalently:
s_mm <- M3X3_seasonal * (1 - M2X12)
#> Error: object 'M3X3_seasonal' not found
s <- y * s_mm
#> Error: object 's_mm' not found
plot(s)
#> Error: object 's' not found
```
