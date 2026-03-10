# Get properties of filters

Get properties of filters

## Usage

``` r
get_properties_function(
  x,
  component = c("Symmetric Gain", "Symmetric Phase", "Symmetric transfer",
    "Asymmetric Gain", "Asymmetric Phase", "Asymmetric transfer"),
  ...
)
```

## Arguments

- x:

  a `"moving_average"` or `"finite_filters"` object.

- component:

  the component to extract.

- ...:

  unused other arguments.

## Examples

``` r
filter <- lp_filter(3, kernel = "Henderson")
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): RcallMethod: cannot determine object class
sgain <- get_properties_function(filter, "Symmetric Gain")
#> Error in UseMethod("get_properties_function", x): no applicable method for 'get_properties_function' applied to an object of class "function"
plot(sgain, xlim= c(0, pi/12))
#> Error: object 'sgain' not found
```
