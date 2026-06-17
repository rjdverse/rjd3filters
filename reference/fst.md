# FST criteria

Compute the Fidelity, Smoothness and Timeliness (FST) criteria

## Usage

``` r
fst(weights, lags, passband = pi/6, ...)
```

## Arguments

- weights:

  either a `"moving_average"` or a numeric vector containing weights.

- lags:

  Lags of the moving average (when `weights` is not a
  `"moving_average"`).

- passband:

  Passband threshold for timeliness criterion.

- ...:

  other unused arguments.

## Value

The values of the 3 criteria, the gain and phase of the associated
filter.

## References

Grun-Rehomme, Michel, Fabien Guggemos, and Dominique Ladiray (2018).
“Asymmetric Moving Averages Minimizing Phase Shift”. In: Handbook on
Seasonal Adjustment,
<https://ec.europa.eu/eurostat/web/products-manuals-and-guidelines/-/ks-gq-18-001>.

## Examples

``` r
filter <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "LC")
fst(filter[, "q=0"])
#>   Fidelity Smoothness Timeliness 
#> 0.38785723 1.27229482 0.03034079 
# To compute the statistics on all filters:
fst(filter)
#>                      q=6          q=5          q=4          q=3          q=2
#> Fidelity    2.038158e-01 1.992509e-01 1.879910e-01 1.810991e-01 0.2011061069
#> Smoothness  8.335318e-03 1.727564e-02 2.094998e-02 9.840559e-03 0.0798744087
#> Timeliness -3.009266e-35 3.379695e-05 8.984811e-05 6.850699e-05 0.0003466228
#>                  q=1        q=0
#> Fidelity   0.2678802 0.38785723
#> Smoothness 0.4331967 1.27229482
#> Timeliness 0.0047965 0.03034079
```
