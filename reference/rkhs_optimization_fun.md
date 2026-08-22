# Optimization Function of Reproducing Kernel Hilbert Space (RKHS) Filters

Export function used to compute the optimal bandwidth of Reproducing
Kernel Hilbert Space (RKHS) filters

## Usage

``` r
rkhs_optimization_fun(
  horizon = 6,
  leads = 0,
  degree = 2,
  kernel = c("Biweight", "Henderson", "Epanechnikov", "Triangular", "Uniform",
    "Triweight"),
  asymmetricCriterion = c("Timeliness", "FrequencyResponse", "Accuracy", "Smoothness"),
  density = c("uniform", "rw"),
  passband = 2 * pi/12
)
```

## Arguments

- horizon:

  horizon (bandwidth) of the symmetric filter.

- leads:

  Leads of the filter (should be positive or 0).

- degree:

  degree of polynomial.

- kernel:

  kernel uses.

- asymmetricCriterion:

  the criteria used to compute the optimal bandwidth. If `"Undefined"`,
  \\m+1\\ is used.

- density:

  hypothesis on the spectral density: `"uniform"` (= white noise, the
  default) or `"rw"` (= random walk).

- passband:

  passband threshold.

## Examples

``` r
graphics::plot(
    rkhs_optimization_fun(
        horizon = 6,
        leads = 0,
        degree = 3,
        asymmetricCriterion = "Timeliness"
    ),
    5.5,
    6 * 3,
    ylab = "Timeliness",
    main = "6X0 filter"
)

graphics::plot(
    rkhs_optimization_fun(
        horizon = 6,
        leads = 1,
        degree = 3,
        asymmetricCriterion = "Timeliness"
    ),
    5.5,
    6 * 3,
    ylab = "Timeliness",
    main = "6X1 filter"
)

graphics::plot(
    rkhs_optimization_fun(
        horizon = 6,
        leads = 2,
        degree = 3,
        asymmetricCriterion = "Timeliness"
    ),
    5.5,
    6 * 3,
    ylab = "Timeliness",
    main = "6X2 filter"
)

graphics::plot(
    rkhs_optimization_fun(
        horizon = 6,
        leads = 3,
        degree = 3,
        asymmetricCriterion = "Timeliness"
    ),
    5.5,
    6 * 3,
    ylab = "Timeliness",
    main = "6X3 filter"
)

graphics::plot(
    rkhs_optimization_fun(
        horizon = 6,
        leads = 4,
        degree = 3,
        asymmetricCriterion = "Timeliness"
    ),
    5.5,
    6 * 3,
    ylab = "Timeliness",
    main = "6X4 filter"
)

graphics::plot(
    rkhs_optimization_fun(
        horizon = 6,
        leads = 5,
        degree = 3,
        asymmetricCriterion = "Timeliness"
    ),
    5.5,
    6 * 3,
    ylab = "Timeliness",
    main = "6X5 filter"
)
```
