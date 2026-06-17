# Get the coefficients of a kernel

Function to get the coefficient associated to a kernel. Those
coefficients are then used to compute the different filters.

## Usage

``` r
get_kernel(
  kernel = c("Henderson", "Uniform", "Triangular", "Epanechnikov", "Parabolic",
    "BiWeight", "TriWeight", "Tricube", "Trapezoidal", "Gaussian"),
  horizon,
  sd_gauss = 0.25
)
```

## Arguments

- kernel:

  kernel uses.

- horizon:

  horizon (bandwidth) of the symmetric filter.

- sd_gauss:

  standard deviation for gaussian kernel. By default 0.25.

## Value

`tskernel` object (see [kernel](https://rdrr.io/r/stats/kernel.html)).

## Examples

``` r
get_kernel("Henderson", horizon = 3)
#> Henderson 
#> coef[-3] = 0.04895
#> coef[-2] = 0.13054
#> coef[-1] = 0.20396
#> coef[ 0] = 0.23310
#> coef[ 1] = 0.20396
#> coef[ 2] = 0.13054
#> coef[ 3] = 0.04895
```
