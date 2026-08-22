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
sgain <- get_properties_function(filter, "Symmetric Gain")
graphics::plot(sgain, xlim= c(0, pi/12))
```
