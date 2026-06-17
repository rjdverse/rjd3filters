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
#> $left
#> $left$ar
#> [1] "1.0000"
#> 
#> $left$sar
#> [1] "1.0000"
#> 
#> $left$diff
#> [1] " -  B + 1.0000"
#> 
#> $left$sdiff
#> [1] " -  B^12 + 1.0000"
#> 
#> 
#> $right
#> $right$ma
#> [1] " - 0.4018 B + 1.0000"
#> 
#> $right$sma
#> [1] " - 0.5569 B^12 + 1.0000"
#> 
#> 
```
