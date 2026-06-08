# Validate the `from`/`to` arguments shared by `transmit()` and `progress()`

Validate the `from`/`to` arguments shared by
[`transmit()`](https://accidda.github.io/hestia/reference/transmit.md)
and
[`progress()`](https://accidda.github.io/hestia/reference/progress.md)

## Usage

``` r
check_to_from(from, to)
```

## Arguments

- from:

  origin compartment, expected to be a single non-empty string.

- to:

  destination compartment(s), expected to be a non-empty character
  vector.
