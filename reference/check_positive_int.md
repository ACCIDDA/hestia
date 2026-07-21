# Validate a positive-integer sampler argument

Used by
[`stan_options()`](https://accidda.github.io/hestia/reference/stan_options.md)
to coerce and check scalar count arguments such as `iter`, `chains`, and
`warmup`.

## Usage

``` r
check_positive_int(x, name)
```

## Arguments

- x:

  the value to validate.

- name:

  the argument name, used in the error message.

## Value

`x` coerced to a length-one integer.
