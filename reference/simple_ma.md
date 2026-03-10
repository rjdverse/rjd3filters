# Simple Moving Average

A simple moving average is a moving average whose coefficients are all
equal and whose sum is 1

## Usage

``` r
simple_ma(order, lags = -trunc((order - 1)/2))
```

## Arguments

- order:

  number of terms of the moving_average

- lags:

  integer indicating the number of lags of the moving average.

## Examples

``` r
# The M2X12 moving average is computed as
(simple_ma(12, -6) + simple_ma(12, -5)) / 2
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
# The M3X3 moving average is computed as
simple_ma(3, -1) ^ 2
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
# The M3X5 moving average is computed as
simple_ma(3, -1) * simple_ma(5, -2)
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
```
