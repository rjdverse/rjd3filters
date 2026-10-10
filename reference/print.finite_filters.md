# Print method for Finite Filters Objects

Printing objects of class `"finite_filters"`.

## Usage

``` r
# S4 method for class 'finite_filters'
print(
  x,
  sfilter = TRUE,
  rfilters = TRUE,
  lfilters = FALSE,
  zero_as_na = FALSE,
  ...
)

# S4 method for class 'finite_filters'
show(object)

# S3 method for class 'finite_filters'
as.matrix(
  x,
  sfilter = TRUE,
  rfilters = TRUE,
  lfilters = FALSE,
  zero_as_na = FALSE,
  ...
)
```

## Arguments

- x, object:

  An object of class `"finite_filters"`.

- sfilter, rfilters, lfilters:

  Logical indicating whether to add the symmetric filter
  (`sfilter = TRUE` by default), the right filters (`rfilters = TRUE` by
  default) or the left filters (`lfilters = FALSE` by default).

- zero_as_na:

  Logical. Asymmetric filters have fewer terms than the symmetric
  filter. If `TRUE`, the missing past/future terms in the matrix
  representation are filled with `NA` instead of `0`. Defaults to
  `FALSE`.

- ...:

  Further arguments to be passed to or from methods.

## Value

The print method returns the argument `x`, invisibly, as for all
[`base::print()`](https://rdrr.io/r/base/print.html) methods.

`as.matrix.finite_filters()` returns a matrix containing the
coefficients of the symmetric and asymmetric filters. The columns are
named according to the number of future points used (for right filters).
For example, column `q=0` corresponds to the asymmetric filter using no
future points, while `q=1` corresponds to the filter using one future
point, and so on. A similar convention with negative values is used for
left filters (when printed).

## Details

The print method just calls the
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) method and prints
the resulting matrix.

## See also

[`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md).

## Examples

``` r
ff_lp <- lp_filter()
as.matrix(ff_lp)
#>             q=6          q=5          q=4          q=3          q=2
#> t-6 -0.01934985 -0.016429821 -0.010992405 -0.008134877 -0.016032761
#> t-5 -0.02786378 -0.025767846 -0.022036255 -0.020190215 -0.024868237
#> t-4  0.00000000  0.001271838  0.003297605  0.004132155  0.002673996
#> t-3  0.06549178  0.065939529  0.066259471  0.066082532  0.067844235
#> t-2  0.14735651  0.146980166  0.145594283  0.144405855  0.149387420
#> t-1  0.21433675  0.213136306  0.210044599  0.207844681  0.216046109
#> t    0.24005716  0.238032623  0.233235092  0.230023684  0.241444975
#> t+1  0.21433675  0.211488120  0.204984764  0.200761868  0.215403021
#> t+2  0.14735651  0.143683794  0.135474614  0.130240227  0.148101243
#> t+3  0.06549178  0.060994971  0.051079966  0.044834091  0.000000000
#> t+4  0.00000000 -0.005320905 -0.016941735  0.000000000  0.000000000
#> t+5 -0.02786378 -0.034008775  0.000000000  0.000000000  0.000000000
#> t+6 -0.01934985  0.000000000  0.000000000  0.000000000  0.000000000
#>              q=1         q=0
#> t-6 -0.042706925 -0.09186038
#> t-5 -0.038631881 -0.05811026
#> t-4  0.001820871  0.01201758
#> t-3  0.079901630  0.11977342
#> t-2  0.174355336  0.24390220
#> t-1  0.253924544  0.35314649
#> t    0.292233930  0.42113096
#> t+1  0.279102495  0.00000000
#> t+2  0.000000000  0.00000000
#> t+3  0.000000000  0.00000000
#> t+4  0.000000000  0.00000000
#> t+5  0.000000000  0.00000000
#> t+6  0.000000000  0.00000000
```
