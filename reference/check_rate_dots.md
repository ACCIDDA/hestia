# Validate the `gamma`-style transition-rate values passed to `progress()`

Each rate value must be either `NA` (fit the parameter) or a positive
numeric (use the supplied value).

## Usage

``` r
check_rate_dots(dots)
```

## Arguments

- dots:

  the unlisted `...` rate arguments from
  [`progress()`](https://accidda.github.io/hestia/reference/progress.md).
