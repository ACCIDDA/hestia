# Validate a single observation-process specification

Each observation is a named numeric vector of detection probabilities,
one entry per compartment, with values in `[0, 1]`.

## Usage

``` r
check_observation(obs, nm = "observation")
```

## Arguments

- obs:

  a single observation specification.

- nm:

  its name, used in error messages.
