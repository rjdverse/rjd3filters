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

  integer specifying how to imput missing date. `nperiod = 1` means
  imputation using last filtered data (1 period backward),
  `nperiod = 12` with monthly data means imputation using last year
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
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
M2X12 <- (simple_ma(12, -6) + simple_ma(12, -5)) / 2
#> Error in jclassName(class, class.loader = class.loader): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/FiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
composite_ma <- M3X3 * M2X12
#> Error: object 'M3X3' not found
# The last 6 points cannot be computed
composite_ma
#> Error: object 'composite_ma' not found
composite_ma * y
#> Error: object 'composite_ma' not found
# they can be computed using the last filtered data
# e.g. to impute the first 3 missing months with last period:
impute_last_obs(composite_ma, n = 3, nperiod = 1) * y
#> Error: object 'composite_ma' not found
# or using the filtered data of the same month in previous year
impute_last_obs(composite_ma, n = 6, nperiod = 12) * y
#> Error: object 'composite_ma' not found
```
