# Validate the `seed` sampler argument passed to [`stan_options()`](https://accidda.github.io/hestia/reference/stan_options.md)

`seed` must be a single value coercible to an integer.

## Usage

``` r
check_seed(seed)
```

## Arguments

- seed:

  the value to validate.

## Value

`seed` coerced to a length-one integer.
