# Package index

## Create Specific Moving Averages

- [`dfa_filter()`](https://rjdverse.github.io/rjd3filters/reference/dfa_filter.md)
  : Direct Filter Approach
- [`fst_filter()`](https://rjdverse.github.io/rjd3filters/reference/fst_filter.md)
  : Estimation of a filter using the Fidelity-Smoothness-Timeliness
  criteria
- [`localpolynomials()`](https://rjdverse.github.io/rjd3filters/reference/localpolynomials.md)
  : Apply Local Polynomials Filters
- [`lp_filter()`](https://rjdverse.github.io/rjd3filters/reference/lp_filter.md)
  : Local Polynomials Filters
- [`simple_ma()`](https://rjdverse.github.io/rjd3filters/reference/simple_ma.md)
  : Simple Moving Average
- [`rkhs_filter()`](https://rjdverse.github.io/rjd3filters/reference/rkhs_filter.md)
  : Reproducing Kernel Hilbert Space (RKHS) Filters
- [`rkhs_kernel()`](https://rjdverse.github.io/rjd3filters/reference/rkhs_kernel.md)
  : Get RKHS kernel function
- [`rkhs_optimal_bw()`](https://rjdverse.github.io/rjd3filters/reference/rkhs_optimal_bw.md)
  : Optimal Bandwidth of Reproducing Kernel Hilbert Space (RKHS) Filters
- [`rkhs_optimization_fun()`](https://rjdverse.github.io/rjd3filters/reference/rkhs_optimization_fun.md)
  : Optimization Function of Reproducing Kernel Hilbert Space (RKHS)
  Filters
- [`get_kernel()`](https://rjdverse.github.io/rjd3filters/reference/get_kernel.md)
  : Get the coefficients of a kernel
- [`mmsre_filter()`](https://rjdverse.github.io/rjd3filters/reference/mmsre_filter.md)
  : Mean Square Revision Error (mmsre) filter

## Stastics on Moving Averages and estimates

- [`cve()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md)
  [`cv()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md)
  [`loocve()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md)
  [`rt()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md)
  [`cp()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md)
  : Diagnostics and goodness of fit of filtered series
- [`diagnostic_matrix()`](https://rjdverse.github.io/rjd3filters/reference/diagnostic_matrix.md)
  : Compute quality criteria for asymmetric filters
- [`fst()`](https://rjdverse.github.io/rjd3filters/reference/fst.md) :
  FST criteria
- [`implicit_forecasts()`](https://rjdverse.github.io/rjd3filters/reference/implicit_forecasts.md)
  : Retrieve implicit forecasts corresponding to the asymmetric filters
- [`underlying_forecasts()`](https://rjdverse.github.io/rjd3filters/reference/underlying_forecasts.md)
  : Retrieve underlying forecasts corresponding to the asymmetric
  filters
- [`get_properties_function()`](https://rjdverse.github.io/rjd3filters/reference/get_properties_function.md)
  : Get properties of filters
- [`mse()`](https://rjdverse.github.io/rjd3filters/reference/mse.md) :
  Accuracy/smoothness/timeliness criteria through spectral decomposition
- [`var_estimator()`](https://rjdverse.github.io/rjd3filters/reference/var_estimator.md)
  : Variance Estimator
- [`df_var()`](https://rjdverse.github.io/rjd3filters/reference/df_var.md)
  : Compute the degrees of freedom for the variance estimator
- [`confint_filter()`](https://rjdverse.github.io/rjd3filters/reference/confint_filter.md)
  : Confidence intervals
- [`polynomial_matrix()`](https://rjdverse.github.io/rjd3filters/reference/polynomial_matrix.md)
  : Create polynomial matrix

## Manipulation of Moving Averages and Finite Filters

- [`filter()`](https://rjdverse.github.io/rjd3filters/reference/filter.md)
  : Linear Filtering on a Time Series
- [`sum(`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `[`( ``*`<moving_average>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `[`( ``*`<moving_average>`*`,`*`<logical>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `[<-`( ``*`<moving_average>`*`,`*`<ANY>`*`,`*`<missing>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`cbind(`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`rbind(`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<moving_average>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<moving_average>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<numeric>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<moving_average>`*`,`*`<missing>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<moving_average>`*`,`*`<missing>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<moving_average>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<moving_average>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<numeric>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<moving_average>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<moving_average>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<numeric>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<ANY>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<moving_average>`*`,`*`<ANY>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `/`( ``*`<moving_average>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `^`( ``*`<moving_average>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<finite_filters>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<moving_average>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<finite_filters>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<ANY>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<finite_filters>`*`,`*`<ANY>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<numeric>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<finite_filters>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<moving_average>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<finite_filters>`*`,`*`<missing>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<finite_filters>`*`,`*`<missing>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<finite_filters>`*`,`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<moving_average>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<finite_filters>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<numeric>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `/`( ``*`<finite_filters>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `^`( ``*`<finite_filters>`*`,`*`<numeric>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `*`( ``*`<finite_filters>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `+`( ``*`<finite_filters>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `-`( ``*`<finite_filters>`*`,`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `[`( ``*`<finite_filters>`*`,`*`<missing>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  [`` `[`( ``*`<finite_filters>`*`,`*`<ANY>`*`)`](https://rjdverse.github.io/rjd3filters/reference/filters_operations.md)
  : Operations on Filters
- [`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md)
  [`is.finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md)
  [`show(`*`<finite_filters>`*`)`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md)
  : Manipulating Finite Filters
- [`get_moving_average()`](https://rjdverse.github.io/rjd3filters/reference/get_moving_average.md)
  : Get Moving Averages from ARIMA model
- [`impute_last_obs()`](https://rjdverse.github.io/rjd3filters/reference/impute_last_obs.md)
  : Impute Incomplete Finite Filters
- [`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`is.moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`is_symmetric()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`upper_bound()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`lower_bound()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`mirror()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`rev(`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`length(`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`to_seasonal()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  [`show(`*`<moving_average>`*`)`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
  : Manipulation of moving averages
- [`plot_coef()`](https://rjdverse.github.io/rjd3filters/reference/plot_filters.md)
  [`plot_gain()`](https://rjdverse.github.io/rjd3filters/reference/plot_filters.md)
  [`plot_phase()`](https://rjdverse.github.io/rjd3filters/reference/plot_filters.md)
  : Plots filters properties

## Data

- [`retailsa`](https://rjdverse.github.io/rjd3filters/reference/retailsa.md)
  : Seasonally Adjusted Retail Sales
