# Impute Incomplete Finite Filters

Impute Incomplete Finite Filters

## Usage

``` r
impute_last_obs(x, n, nperiod = 1, backward = TRUE, forward = TRUE)
```

## Arguments

- x:

  a
  [`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md)
  object.

- n:

  integer specifying the number of imputed periods. By default all the
  missing moving averages are imputed.

- nperiod:

  integer specifying how to impute missing date.

  - `nperiod = 1` means imputation using last filtered data (1 period
    backward),

  - `nperiod = 12` with monthly data means imputation using last year
    filtered data, etc.

- backward, forward:

  boolean indicating if the imputation should be done backward (on left
  filters), forward (on right filters).

## Details

When combining finite filters and a moving average, the first and/or the
last points cannot be computed.

For example, using the M2X12 moving average, that is to say the
symmetric moving average with coefficients \$\$ \theta =
\frac{1}{24}B^{6} +
\frac{1}{12}B^{5}+\dots+\frac{1}{12}B^{-5}+\frac{1}{24}B^{-6}, \$\$ the
first and last 6 points cannot be computed.

`impute_last_obs()` allows to impute the first/last points using the
`nperiod` previous filtered data. With `nperiod = 1`, the last filtered
data is used for the imputation, with `nperiod = 12` and monthly data,
the last year filtered data is used for the imputation, etc.

## Examples

``` r
y <- window(retailsa$AllOtherGenMerchandiseStores, start = 2008)
M3 <- moving_average(rep(1/3, 3), lags = -1)
M3X3 <- M3 * M3
M2X12 <- (simple_ma(12, -6) + simple_ma(12, -5)) / 2
composite_ma <- M3X3 * M2X12
# The last 6 points cannot be computed
composite_ma
#> [1] "0.0046 B^8 + 0.0185 B^7 + 0.0417 B^6 + 0.0648 B^5 + 0.0787 B^4 + 0.0833 B^3 + 0.0833 B^2 + 0.0833 B + 0.0833 + 0.0833 F + 0.0833 F^2 + 0.0833 F^3 + 0.0787 F^4 + 0.0648 F^5 + 0.0417 F^6 + 0.0185 F^7 + 0.0046 F^8"
composite_ma * y
#>           Jan      Feb      Mar      Apr      May      Jun      Jul      Aug
#> 2008       NA       NA       NA       NA       NA       NA       NA       NA
#> 2009 3874.468 3887.489 3906.257 3930.825 3959.576 3991.101 4025.292 4062.564
#> 2010 4253.845 4289.213 4322.790 4353.753       NA       NA       NA       NA
#>           Sep      Oct      Nov      Dec
#> 2008 3835.187 3845.520 3855.656 3864.796
#> 2009 4102.100 4142.065 4180.559 4217.600
#> 2010       NA       NA       NA       NA
# they can be computed using the last filtered data
# e.g. to impute the first 3 missing months with last period:
impute_last_obs(composite_ma, n = 3, nperiod = 1) * y
#>           Jan      Feb      Mar      Apr      May      Jun      Jul      Aug
#> 2008       NA       NA       NA       NA       NA 3835.187 3835.187 3835.187
#> 2009 3874.468 3887.489 3906.257 3930.825 3959.576 3991.101 4025.292 4062.564
#> 2010 4253.845 4289.213 4322.790 4353.753 4353.753 4353.753 4353.753       NA
#>           Sep      Oct      Nov      Dec
#> 2008 3835.187 3845.520 3855.656 3864.796
#> 2009 4102.100 4142.065 4180.559 4217.600
#> 2010       NA       NA       NA       NA
# or using the filtered data of the same month in previous year
impute_last_obs(composite_ma, n = 6, nperiod = 12) * y
#>           Jan      Feb      Mar      Apr      May      Jun      Jul      Aug
#> 2008       NA       NA 3906.257 3930.825 3959.576 3991.101 4025.292 4062.564
#> 2009 3874.468 3887.489 3906.257 3930.825 3959.576 3991.101 4025.292 4062.564
#> 2010 4253.845 4289.213 4322.790 4353.753 3959.576 3991.101 4025.292 4062.564
#>           Sep      Oct      Nov      Dec
#> 2008 3835.187 3845.520 3855.656 3864.796
#> 2009 4102.100 4142.065 4180.559 4217.600
#> 2010 4102.100 4142.065       NA       NA
```
