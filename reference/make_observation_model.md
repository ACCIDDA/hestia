# Make Observation Process Model

Make Observation Process Model

## Usage

``` r
make_observation_model(...)
```

## Arguments

- ...:

  A series of named vectors. Each vector corresponds to an observation
  type. Each entry in the vector is named for the corresponding
  compartment in the infection process model. The value is the
  probability of observing a positive observation given the individual
  is in the compartment.

## Value

A list with one entry for each observation type.

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
```
