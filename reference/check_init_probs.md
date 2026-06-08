# Validate the `init_probs` argument passed to `run_model()`

Initial state probabilities must be a finite non-negative numeric vector
that sums to roughly one.

## Usage

``` r
check_init_probs(init_probs)
```

## Arguments

- init_probs:

  the vector of initial state probabilities.
