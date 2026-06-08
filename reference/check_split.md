# Validate the `split` argument shared by `transmit()` and `progress()`

`split` may be left as a single `NA` (no split), a character vector of
parameter names, or a numeric vector of proportions. This only checks
the broad type and finiteness; the detailed length/sum rules stay in
[`split_check()`](https://accidda.github.io/hestia/reference/split_check.md).

## Usage

``` r
check_split(split)
```

## Arguments

- split:

  the split specification to validate.
