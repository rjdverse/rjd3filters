# Accuracy/smoothness/timeliness criteria through spectral decomposition

Accuracy/smoothness/timeliness criteria through spectral decomposition

## Usage

``` r
mse(aweights, sweights, density = c("uniform", "rw"), passband = pi/6, ...)
```

## Arguments

- aweights:

  `moving_average` object or weights of the asymmetric filter (from -n
  to m).

- sweights:

  `moving_average` object or weights of the symmetric filter (from 0 to
  n or -n to n).

- density:

  hypothesis on the spectral density: `"uniform"` (= white noise, the
  default) or `"rw"` (= random walk).

- passband:

  passband threshold.

- ...:

  other unused arguments.

## Value

The criteria

## References

Wildi, Marc and McElroy, Tucker (2019). “The trilemma between accuracy,
timeliness and smoothness in real-time signal extraction”. In:
International Journal of Forecasting 35.3, pp. 1072–1084.

## Examples

``` r
filter <- lp_filter(horizon = 6, kernel = "Henderson", endpoints = "LC")
sweights <- filter[, "q=6"]
aweights <- filter[, "q=0"]
mse(aweights, sweights)
#>   accuracy smoothness timeliness   residual 
#> 0.01507927 0.52517037 0.05226739 0.31059437 
# Or to compute directly the criteria on all asymmetric filters:
mse(filter)
#>                      q=6          q=5          q=4          q=3          q=2
#> accuracy    6.024615e-31 3.049917e-04 0.0023586820 0.0039385492 0.0010541330
#> smoothness  2.495366e-32 1.407955e-03 0.0051133637 0.0059461839 0.0181732313
#> timeliness  3.783119e-31 6.949349e-05 0.0001954977 0.0001533431 0.0007302353
#> residual   -4.091914e-31 1.376327e-03 0.0052601574 0.0050943353 0.0144950219
#>                    q=1        q=0
#> accuracy   0.001615599 0.01507927
#> smoothness 0.147090323 0.52517037
#> timeliness 0.009039553 0.05226739
#> residual   0.076237423 0.31059437
```
