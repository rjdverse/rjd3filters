# Operations on Filters

Manipulation of
[`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
or
[`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md)
objects.

## Usage

``` r
# S3 method for class 'moving_average'
sum(..., na.rm = FALSE)

# S4 method for class 'moving_average,numeric'
x[i]

# S4 method for class 'moving_average,logical'
x[i]

# S4 method for class 'moving_average,ANY,missing,numeric'
x[i] <- value

# S3 method for class 'moving_average'
cbind(..., zero_as_na = FALSE)

# S3 method for class 'moving_average'
rbind(...)

# S4 method for class 'moving_average,moving_average'
e1 + e2

# S4 method for class 'moving_average,numeric'
e1 + e2

# S4 method for class 'numeric,moving_average'
e1 + e2

# S4 method for class 'moving_average,missing'
e1 + e2

# S4 method for class 'moving_average,missing'
e1 - e2

# S4 method for class 'moving_average,moving_average'
e1 - e2

# S4 method for class 'moving_average,numeric'
e1 - e2

# S4 method for class 'numeric,moving_average'
e1 - e2

# S4 method for class 'moving_average,moving_average'
e1 * e2

# S4 method for class 'moving_average,numeric'
e1 * e2

# S4 method for class 'numeric,moving_average'
e1 * e2

# S4 method for class 'ANY,moving_average'
e1 * e2

# S4 method for class 'moving_average,ANY'
e1 * e2

# S4 method for class 'moving_average,numeric'
e1/e2

# S4 method for class 'moving_average,numeric'
e1^e2

# S4 method for class 'finite_filters,moving_average'
e1 * e2

# S4 method for class 'moving_average,finite_filters'
e1 * e2

# S4 method for class 'finite_filters,numeric'
e1 * e2

# S4 method for class 'ANY,finite_filters'
e1 * e2

# S4 method for class 'finite_filters,ANY'
e1 * e2

# S4 method for class 'numeric,finite_filters'
e1 + e2

# S4 method for class 'finite_filters,moving_average'
e1 + e2

# S4 method for class 'moving_average,finite_filters'
e1 + e2

# S4 method for class 'finite_filters,missing'
e1 + e2

# S4 method for class 'finite_filters,missing'
e1 - e2

# S4 method for class 'finite_filters,moving_average'
e1 - e2

# S4 method for class 'moving_average,finite_filters'
e1 - e2

# S4 method for class 'finite_filters,numeric'
e1 - e2

# S4 method for class 'numeric,finite_filters'
e1 - e2

# S4 method for class 'finite_filters,numeric'
e1/e2

# S4 method for class 'finite_filters,numeric'
e1^e2

# S4 method for class 'finite_filters,finite_filters'
e1 * e2

# S4 method for class 'finite_filters,finite_filters'
e1 + e2

# S4 method for class 'finite_filters,finite_filters'
e1 - e2

# S4 method for class 'finite_filters,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'finite_filters,ANY'
x[i, j, ..., drop = TRUE]
```

## Arguments

- ..., drop, na.rm:

  other parameters.

- x, e1, e2:

  object.

- i, j, value:

  indices specifying elements to extract or replace and the new value

- zero_as_na:

  boolean indicating if, when merging several moving averages (`cbind`)
  if trailing and leading zeros added to have a matrix form should be
  replaced by `NA`.

## Value

A
[`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
or
[`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md)
object.

## Examples

``` r
e1 <- moving_average(rep(1,12), lags = -6)
e1 <- e1/sum(e1)
e2 <- moving_average(rep(1/12, 12), lags = -5)
M2X12 <- (e1 + e2)/2
M2X12
#> [1] "0.0417 B^6 + 0.0833 B^5 + 0.0833 B^4 + 0.0833 B^3 + 0.0833 B^2 + 0.0833 B + 0.0833 + 0.0833 F + 0.0833 F^2 + 0.0833 F^3 + 0.0833 F^4 + 0.0833 F^5 + 0.0417 F^6"
h13 <- lp_filter(horizon = 6)
h13
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
h13 * M2X12
#>                q=6           q=5           q=4           q=3           q=2
#> t-12 -0.0008062436 -0.0006845759 -0.0004580169 -0.0003389532 -0.0006680317
#> t-11 -0.0027734778 -0.0024428120 -0.0018342111 -0.0015191654 -0.0023722399
#> t-10 -0.0039344685 -0.0034634790 -0.0026149881 -0.0021882512 -0.0032970000
#> t-9  -0.0012056442 -0.0006630054  0.0002832234  0.0007373607 -0.0003587404
#> t-8   0.0076630348  0.0082086486  0.0091104631  0.0095077102  0.0086925786
#> t-7   0.0227335874  0.0232135016  0.0239287499  0.0241848158  0.0239189756
#> t-6   0.0416666667  0.0420122070  0.0423987370  0.0424293310  0.0429811041
#> t-5   0.0605997460  0.0607422379  0.0606578977  0.0603787290  0.0620164372
#> t-4   0.0756702985  0.0755410677  0.0748437051  0.0741704830  0.0771624482
#> t-3   0.0845389775  0.0840693496  0.0826168126  0.0814652462  0.0833333333
#> t-2   0.0872678019  0.0863891023  0.0840392389  0.0833333333  0.0833333333
#> t-1   0.0861068111  0.0847503656  0.0833333333  0.0833333333  0.0833333333
#> t     0.0849458204  0.0840179092  0.0837913502  0.0836722865  0.0840013650
#> t+1   0.0861068111  0.0857761453  0.0851675444  0.0848524987  0.0857055733
#> t+2   0.0872678019  0.0867968123  0.0859483215  0.0855215846  0.0866303333
#> t+3   0.0845389775  0.0839963387  0.0830501100  0.0825959726  0.0836920737
#> t+4   0.0756702985  0.0751246847  0.0742228702  0.0738256232  0.0746407548
#> t+5   0.0605997460  0.0601198318  0.0594045834  0.0591485175  0.0594143577
#> t+6   0.0416666667  0.0413211264  0.0409345963  0.0409040023  0.0403522292
#> t+7   0.0227335874  0.0225910954  0.0226754356  0.0229546043  0.0213168961
#> t+8   0.0076630348  0.0077922656  0.0084896282  0.0091628504  0.0061708851
#> t+9  -0.0012056442 -0.0007360163  0.0007165207  0.0018680871  0.0000000000
#> t+10 -0.0039344685 -0.0030557690 -0.0007059056  0.0000000000  0.0000000000
#> t+11 -0.0027734778 -0.0014170323  0.0000000000  0.0000000000  0.0000000000
#> t+12 -0.0008062436  0.0000000000  0.0000000000  0.0000000000  0.0000000000
#>               q=1          q=0
#> t-12 -0.001779455 -0.003827516
#> t-11 -0.005168572 -0.010076292
#> t-10 -0.006702364 -0.011996821
#> t-9  -0.003297260 -0.006505530
#> t-8   0.007296780  0.008647621
#> t-7   0.025141775  0.033524650
#> t-6   0.047898378  0.065786210
#> t-5   0.071704063  0.083333333
#> t-4   0.083333333  0.083333333
#> t-3   0.083333333  0.083333333
#> t-2   0.083333333  0.083333333
#> t-1   0.083333333  0.083333333
#> t     0.085112789  0.087160849
#> t+1   0.088501905  0.093409626
#> t+2   0.090035698  0.095330154
#> t+3   0.086630593  0.089838863
#> t+4   0.076036553  0.074685712
#> t+5   0.058191558  0.049808683
#> t+6   0.035434955  0.017547123
#> t+7   0.011629271  0.000000000
#> t+8   0.000000000  0.000000000
#> t+9   0.000000000  0.000000000
#> t+10  0.000000000  0.000000000
#> t+11  0.000000000  0.000000000
#> t+12  0.000000000  0.000000000
```
