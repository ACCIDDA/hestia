# Estimate Test Sensitivity and Specificity

Specifies an observation type whose sensitivity and specificity are
estimated rather than fixed, for use in
[make_observation_model](https://accidda.github.io/hestia/reference/make_observation_model.md).
Compartments in `positive` give a positive result with probability equal
to the sensitivity, and all other compartments with probability one
minus the specificity. Sensitivity is constrained to exceed one minus
specificity, which rules out the mirror-image solution in which the test
is read the wrong way round.

## Usage

``` r
estimate_test(positive, se, sp)
```

## Arguments

- positive:

  Character vector of compartments in which the target is present, so
  the test should be positive.

- se:

  Numeric vector `c(a, b)` giving the Beta(a, b) prior on sensitivity.

- sp:

  Numeric vector `c(a, b)` giving the Beta(a, b) prior on specificity.

## Value

A list of class `hestia_estimated_test`.

## Examples

``` r
estimate_test(positive = "I", se = c(9, 1), sp = c(19, 1))
#> $positive
#> [1] "I"
#> 
#> $se
#> [1] 9 1
#> 
#> $sp
#> [1] 19  1
#> 
#> attr(,"class")
#> [1] "hestia_estimated_test"
```
