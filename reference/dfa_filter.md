# Direct Filter Approach

Direct Filter Approach

## Usage

``` r
dfa_filter(
  horizon = 6,
  degree = 0,
  density = c("uniform", "rw"),
  targetfilter = lp_filter(horizon = horizon)@sfilter,
  passband = 2 * pi/12,
  accuracy.weight = 1/3,
  smoothness.weight = 1/3,
  timeliness.weight = 1/3
)
```

## Arguments

- horizon:

  horizon (bandwidth) of the symmetric filter.

- degree:

  degree of polynomial.

- density:

  hypothesis on the spectral density: `"uniform"` (= white noise, the
  default) or `"rw"` (= random walk).

- targetfilter:

  the weights of the symmetric target filters (by default the Henderson
  filter).

- passband:

  passband threshold.

- accuracy.weight, smoothness.weight, timeliness.weight:

  the weight used for the optimisation. The weight associated to the
  residual is derived so that the sum of the four weights are equal to
  1.

## Details

Moving average computed by a minimisation of a weighted mean of three
criteria under polynomials constraints. The criteria come from the
decomposition of the mean squared error between th trend-cycle

Let \\\boldsymbol \theta=(\theta\_{-p},\dots,\theta\_{f})'\\ be a moving
average where \\p\\ and \\f\\ are two integers defined by the parameter
`lags` and `leads`. The three criteria are:

## Examples

``` r
# \donttest{
dfa_filter(horizon = 6, degree = 0)
#>             q=6          q=5           q=4         q=3         q=2         q=1
#> t-6 -0.01934985 -0.030012305 -0.0408092039 -0.04068978 -0.03622977 -0.03672904
#> t-5 -0.02786378 -0.033359743 -0.0411977083 -0.04106306 -0.03072500 -0.01681858
#> t-4  0.00000000  0.001617955  0.0031630124  0.00320233  0.01204228  0.04093710
#> t-3  0.06549178  0.071870402  0.0795170526  0.07941604  0.08550769  0.12220301
#> t-2  0.14735651  0.154367848  0.1632402190  0.16312197  0.16726498  0.20144314
#> t-1  0.21433675  0.218278871  0.2226926829  0.22264487  0.22680122  0.25014591
#> t    0.24005716  0.238829368  0.2358094694  0.23588451  0.24196319  0.24750757
#> t+1  0.21433675  0.207944469  0.1985282507  0.19868118  0.20515234  0.19131090
#> t+2  0.14735651  0.138054450  0.1272631817  0.12741087  0.12822306  0.00000000
#> t+3  0.06549178  0.057433857  0.0513680901  0.05139107  0.00000000  0.00000000
#> t+4  0.00000000 -0.002016302  0.0004249535  0.00000000  0.00000000  0.00000000
#> t+5 -0.02786378 -0.023008869  0.0000000000  0.00000000  0.00000000  0.00000000
#> t+6 -0.01934985  0.000000000  0.0000000000  0.00000000  0.00000000  0.00000000
#>             q=0
#> t-6 -0.02675707
#> t-5  0.01420277
#> t-4  0.08880802
#> t-3  0.17595533
#> t-2  0.24508791
#> t-1  0.26838604
#> t    0.23431701
#> t+1  0.00000000
#> t+2  0.00000000
#> t+3  0.00000000
#> t+4  0.00000000
#> t+5  0.00000000
#> t+6  0.00000000
dfa_filter(horizon = 6, degree = 2)
#>             q=6          q=5         q=4         q=3         q=2          q=1
#> t-6 -0.01934985 -0.031105441 -0.04888338 -0.05764640 -0.04948784 -0.007961043
#> t-5 -0.02786378 -0.030840929 -0.03722806 -0.04096008 -0.03686433 -0.036884447
#> t-4  0.00000000  0.005553125  0.01509948  0.02043095  0.01419427 -0.019710732
#> t-3  0.06549178  0.075680118  0.09443554  0.10664294  0.09370622  0.043998864
#> t-2  0.14735651  0.157683962  0.17753428  0.18829824  0.17738803  0.138819404
#> t-1  0.21433675  0.221376079  0.23306525  0.23716685  0.23398281  0.237104542
#> t    0.24005716  0.241802252  0.23989975  0.23478728  0.24355756  0.308924895
#> t+1  0.21433675  0.210178737  0.19487115  0.18311315  0.20164062  0.335708518
#> t+2  0.14735651  0.138186344  0.11635204  0.10287717  0.12188266  0.000000000
#> t+3  0.06549178  0.054026865  0.03459332  0.02528989  0.00000000  0.000000000
#> t+4  0.00000000 -0.009233115 -0.01973936  0.00000000  0.00000000  0.000000000
#> t+5 -0.02786378 -0.033307998  0.00000000  0.00000000  0.00000000  0.000000000
#> t+6 -0.01934985  0.000000000  0.00000000  0.00000000  0.00000000  0.000000000
#>             q=0
#> t-6  0.13737578
#> t-5 -0.08485848
#> t-4 -0.16654786
#> t-3 -0.10200794
#> t-2  0.09325916
#> t-1  0.38573469
#> t    0.73704465
#> t+1  0.00000000
#> t+2  0.00000000
#> t+3  0.00000000
#> t+4  0.00000000
#> t+5  0.00000000
#> t+6  0.00000000
# }
```
