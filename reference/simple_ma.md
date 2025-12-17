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
#> [1] "0.0417 B^6 + 0.0833 B^5 + 0.0833 B^4 + 0.0833 B^3 + 0.0833 B^2 + 0.0833 B + 0.0833 + 0.0833 F + 0.0833 F^2 + 0.0833 F^3 + 0.0833 F^4 + 0.0833 F^5 + 0.0417 F^6"
# The M3X3 moving average is computed as
simple_ma(3, -1) ^ 2
#> [1] "0.1111 B^2 + 0.2222 B + 0.3333 + 0.2222 F + 0.1111 F^2"
# The M3X5 moving average is computed as
simple_ma(3, -1) * simple_ma(5, -2)
#> [1] "0.0667 B^3 + 0.1333 B^2 + 0.2000 B + 0.2000 + 0.2000 F + 0.1333 F^2 + 0.0667 F^3"
```
