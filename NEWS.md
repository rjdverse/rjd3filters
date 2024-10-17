# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Added

* New function `polynomial_matrix()` to create a matrix of polynomial terms.

* New function `mmsre_filter()` to compute the general Proietti and Luati (2008) filter with extension for non symmetric filters and with Timeliness criterion.

### Changed

* `filter()` correction when the length of the series equals the length of the filter.


## [2.1.1] - 2024-12-07

### Changed

* `cross_validation()` function renamed to `cve()`, `cross_validation()` is now deprecated.
* New .jar files added (related to [v2.2.0](https://github.com/jdemetra/jdplus-incubator/releases/tag/v2.2.0))

### Added

* New functions to compute functions to compute diagnostics and goodness of fit of filtered series: cross validation (`cv()`) and cross validate estimate (`cve()`), leave-one-out cross validation estimate (`loocve`), CP statistic (`cp()`) and Rice's T statistics (`rt()`).

* New function `confint_filter()` to compute confidence intervals for filtered series.

* New function `is.finite_filters()`.

* New parameter `zero_as_na` in `cbind.moving_average`, boolean indicating if trealing and leading zeros added to have a matrix form should be replaced by `NA`.


## [2.1.0] - 2024-04-18

### Changed

* Merge pull request #22 from rjdemetra/develop
* v2.1.0 - last improvements


## [2.0.0] - 2023-12-12

### Changed

* Merge pull request #12 from rjdemetra/develop
* Merge pull request #11 from rjdemetra/main


## [1.0.0] - 2023-07-06

### Added

* New Jars


[Unreleased]: https://github.com/rjdverse/rjd3filters/compare/v2.1.1...HEAD
[2.1.1]: https://github.com/rjdverse/rjd3filters/releases/tag/v2.1.0...v2.1.1
[2.1.0]: https://github.com/rjdverse/rjd3filters/releases/tag/v2.0.0...v2.1.0
[2.0.0]: https://github.com/rjdverse/rjd3filters/releases/tag/v1.0.0...v2.0.0
[1.0.0]: https://github.com/rjdverse/rjd3filters/releases/tag/v1.0.0
