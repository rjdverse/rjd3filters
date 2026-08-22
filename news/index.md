# Changelog

## rjd3filters 2.4.0.9000

All notable changes to this project will be documented in this file.

The format is based on [Keep a
Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### [Unreleased](https://github.com/rjdverse/rjd3filters/compare/v2.3.0...HEAD)

#### Added

- New method
  [`underlying_forecasts()`](https://rjdverse.github.io/rjd3filters/reference/underlying_forecasts.md).
- New JARS related to version
  [3.7.1](https://github.com/jdemetra/jdplus-main/releases/tag/v3.7.1)
  of JDemetra+.

#### Changed

- Correction in internal functions when filtering series with NA and the
  beginning and not and the end of the series.
- For polynomial methods, default I/C ratio fixed to 3.5 (as in X-11 for
  H-13).
- `implicit_forecast()` function renamed to
  [`implicit_forecasts()`](https://rjdverse.github.io/rjd3filters/reference/implicit_forecasts.md).

### [2.3.0](https://github.com/rjdverse/rjd3filters/compare/v2.2.0...v2.3.0) - 2025-04-24

#### Changed

- New JARS related to version
  [2.3.0](https://github.com/jdemetra/jdplus-incubator/releases/tag/v2.3.0)

### [2.2.0](https://github.com/rjdverse/rjd3filters/compare/v2.1.1...v2.2.0) - 2025-03-01

#### Added

- New function
  [`polynomial_matrix()`](https://rjdverse.github.io/rjd3filters/reference/polynomial_matrix.md)
  to create a matrix of polynomial terms.
- New function
  [`mmsre_filter()`](https://rjdverse.github.io/rjd3filters/reference/mmsre_filter.md)
  to compute the general Proietti and Luati (2008) filter with extension
  for non symmetric filters and with Timeliness criterion.
- New parameter to
  [`confint_filter()`](https://rjdverse.github.io/rjd3filters/reference/confint_filter.md)
  to specify if the variance should be estimated for each asymmetric
  filters (default) instead of using the variance associated the
  symmetric estimates.

#### Changed

- [`filter()`](https://rjdverse.github.io/rjd3filters/reference/filter.md)
  correction when the length of the series equals the length of the
  filter.
- [`confint_filter()`](https://rjdverse.github.io/rjd3filters/reference/confint_filter.md)
  uses by default a Student distribution instead of a Normal
  distribution.

### [2.1.1](https://github.com/rjdverse/rjd3filters/compare/v2.1.0...v2.1.1) - 2024-07-12

#### Added

- New functions to compute functions to compute diagnostics and goodness
  of fit of filtered series: cross validation
  ([`cv()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md))
  and cross validate estimate
  ([`cve()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md)),
  leave-one-out cross validation estimate (`loocve`), CP statistic
  ([`cp()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md))
  and Rice’s T statistics
  ([`rt()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md)).
- New function
  [`confint_filter()`](https://rjdverse.github.io/rjd3filters/reference/confint_filter.md)
  to compute confidence intervals for filtered series.
- New function
  [`is.finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md).
- New parameter `zero_as_na` in `cbind.moving_average`, boolean
  indicating if trealing and leading zeros added to have a matrix form
  should be replaced by `NA`.

#### Changed

- `cross_validation()` function renamed to
  [`cve()`](https://rjdverse.github.io/rjd3filters/reference/diagnostics-fit.md),
  `cross_validation()` is now deprecated.
- New JARS related to version
  [2.2.0](https://github.com/jdemetra/jdplus-incubator/releases/tag/v2.2.0)

### [2.1.0](https://github.com/rjdverse/rjd3filters/compare/v2.0.0...v2.1.0) - 2024-04-18

#### Changed

- Merge pull request
  [\#22](https://github.com/rjdverse/rjd3filters/issues/22) from
  rjdemetra/develop
- New JARS related to version
  [2.1.0](https://github.com/jdemetra/jdplus-incubator/releases/tag/v2.1.0)

### [2.0.0](https://github.com/rjdverse/rjd3filters/compare/v1.0.0...v2.0.0) - 2023-12-12

#### Changed

- Merge pull request
  [\#12](https://github.com/rjdverse/rjd3filters/issues/12) from
  rjdemetra/develop
- Merge pull request
  [\#11](https://github.com/rjdverse/rjd3filters/issues/11) from
  rjdemetra/main
- New JARS related to version
  [2.0.0](https://github.com/jdemetra/jdplus-incubator/releases/tag/v2.0.0)

### [1.0.0](https://github.com/rjdverse/rjd3filters/releases/tag/v1.0.0) - 2023-07-06

#### Added

- Initial JARS related to version
  [1.0.0](https://github.com/jdemetra/jdplus-incubator/releases/tag/v1.0.0)
