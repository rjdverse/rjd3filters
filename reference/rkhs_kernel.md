# Get RKHS kernel function

Get RKHS kernel function

## Usage

``` r
rkhs_kernel(
  kernel = c("Biweight", "Henderson", "Epanechnikov", "Triangular", "Uniform",
    "Triweight"),
  degree = 2,
  horizon = 6
)
```

## Arguments

- kernel:

  kernel uses.

- degree:

  degree of polynomial.

- horizon:

  horizon (bandwidth) of the symmetric filter.

## Value

A function that takes a numeric input and returns the value of the RKHS
kernel.

## Examples

``` r
biweight <- rkhs_kernel(kernel = "Biweight")
triangular <- rkhs_kernel(kernel = "Triangular")
graphics::plot(biweight, -1, 1)
graphics::plot(triangular, -1, 1, add = TRUE, col = "orange")
```
