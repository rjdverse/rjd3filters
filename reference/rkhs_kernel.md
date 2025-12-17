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
