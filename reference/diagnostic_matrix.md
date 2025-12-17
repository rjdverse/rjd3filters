# Compute quality criteria for asymmetric filters

Function du compute a diagnostic matrix of quality criteria for
asymmetric filters

## Usage

``` r
diagnostic_matrix(x, lags, passband = pi/6, sweights, ...)
```

## Arguments

- x:

  Weights of the asymmetric filter (from -lags to m).

- lags:

  Lags of the filter (should be positive).

- passband:

  passband threshold.

- sweights:

  Weights of the symmetric filter (from 0 to lags or -lags to lags). If
  missing, the criteria from the functions
  [`mse`](https://rjdverse.github.io/rjd3filters/reference/mse.md) are
  not computed.

- ...:

  optional arguments to
  [`mse`](https://rjdverse.github.io/rjd3filters/reference/mse.md).

## Details

For a moving average of coefficients \\\theta=(\theta_i)\_{-p\le i\le
q}\\ `diagnostic_matrix` returns a `list` with the following ten
criteria:

- `b_c` Constant bias (if \\b_c=0\\, \\\theta\\ preserve constant
  trends) \$\$\sum\_{i=-p}^q\theta_i - 1\$\$

- `b_l` Linear bias (if \\b_c=b_l=0\\, \\\theta\\ preserve constant
  trends) \$\$\sum\_{i=-p}^q i \theta_i\$\$

- `b_q` Quadratic bias (if \\b_c=b_l=b_q=0\\, \\\theta\\ preserve
  quadratic trends) \$\$\sum\_{i=-p}^q i^2 \theta_i\$\$

- `F_g` Fidelity criterium of Grun-Rehomme et al (2018) \$\$\$\$

- `S_g` Smoothness criterium of Grun-Rehomme et al (2018)

- `T_g` Timeliness criterium of Grun-Rehomme et al (2018)

- `A_w` Accuracy criterium of Wildi and McElroy (2019)

- `S_w` Smoothness criterium of Wildi and McElroy (2019)

- `T_w` Timeliness criterium of Wildi and McElroy (2019)

- `R_w` Residual criterium of Wildi and McElroy (2019)

## References

Grun-Rehomme, Michel, Fabien Guggemos, and Dominique Ladiray (2018).
“Asymmetric Moving Averages Minimizing Phase Shift”. In: Handbook on
Seasonal Adjustment.

Wildi, Marc and McElroy, Tucker (2019). “The trilemma between accuracy,
timeliness and smoothness in real-time signal extraction”. In:
International Journal of Forecasting 35.3, pp. 1072–1084.
