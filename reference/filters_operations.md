# Operations on Filters

Manipulation of
[`moving_average()`](https://rjdverse.github.io/rjd3filters/reference/moving_average.md)
or
[`finite_filters()`](https://rjdverse.github.io/rjd3filters/reference/finite_filters.md)
objects

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

  object

- i, j, value:

  indices specifying elements to extract or replace and the new value

- zero_as_na:

  boolean indicating if, when merging several moving averages (`cbind`)
  if trailing and leading zeros added to have a matrix form should be
  replaced by `NA`.
