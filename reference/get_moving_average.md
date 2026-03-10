# Get Moving Averages from ARIMA model

Get Moving Averages from ARIMA model

## Usage

``` r
get_moving_average(x, ...)
```

## Arguments

- x:

  the object.

- ...:

  unused parameters

## Examples

``` r
fit <- stats::arima(log10(AirPassengers), c(0, 1, 1),
seasonal = list(order = c(0, 1, 1), period = 12))
get_moving_average(fit)
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
```
