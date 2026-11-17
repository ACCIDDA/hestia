# Make Observation Process Model

Specifies, for each observation type (test), the probability of a
positive result in each compartment. A test's characteristics are either
fixed, or its sensitivity and specificity are estimated with
[estimate_test](https://accidda.github.io/hestia/reference/estimate_test.md).

## Usage

``` r
make_observation_model(...)
```

## Arguments

- ...:

  A series of named observation types. Each is either a named numeric
  vector or an
  [estimate_test](https://accidda.github.io/hestia/reference/estimate_test.md)
  specification. In a numeric vector, each entry is named for the
  corresponding compartment in the infection process model and gives the
  fixed probability of a positive observation given the individual is in
  that compartment.

## Value

A list with one entry for each observation type: a matrix of negative
and positive probabilities for a fixed test, or the
[estimate_test](https://accidda.github.io/hestia/reference/estimate_test.md)
specification for an estimated one.

## Examples

``` r
# Observation process for an SIR model with two observations, PCR test
# results and presence of IgG antibodies.
make_observation_model(
  pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
  igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
)
#> $pcr
#>            S    I    R
#> neg_obs 0.95 0.05 0.95
#> pos_obs 0.05 0.95 0.05
#> 
#> $igg
#>            S    I   R
#> neg_obs 0.99 0.99 0.2
#> pos_obs 0.01 0.01 0.8
#> 

# Estimate PCR sensitivity and specificity, keeping IgG fixed.
make_observation_model(
  pcr = estimate_test(positive = "I", se = c(9, 1), sp = c(19, 1)),
  igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
)
#> $pcr
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
#> 
#> $igg
#>            S    I   R
#> neg_obs 0.99 0.99 0.2
#> pos_obs 0.01 0.01 0.8
#> 
```
