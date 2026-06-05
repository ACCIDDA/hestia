# Build Infection Process Model

Build Infection Process Model

## Usage

``` r
make_infection_model(..., mult_inf_probs = FALSE)
```

## Arguments

- ...:

  a series of
  [progress](https://accidda.github.io/hestia/reference/progress.md) or
  [transmit](https://accidda.github.io/hestia/reference/transmit.md)
  function calls

- mult_inf_probs:

  If FALSE then all infection probabilities are shared across infectious
  compartments

## Value

A data frame with a row for each unique transition in the infection
process model and the following columns:

- "from", the name of the origin compartment

- "to", the name of the destination compartment

- "source", NA (so output columns match with `transmit`)

- "rate_name", name of transition rate parameter if fitting

- "rate_value", numeric value of transition rate if not fitting

- "split_name", name of split parameter if fitting

- "split_value", numeric value of split if not fitting

- "source", names of compartment(s) that are the source of new
  infections

- "mult_inf_probs", a logical indicator for whether infectious
  compartments have separate (`true`) or shared (`FALSE`)
  intra-household infection probabilities

## Examples

``` r
# Basic SIR model with recovery rate gamma to be fit
make_infection_model(transmit(from = "S", to = "I"),
                     progress(from = "I", to = "R", gamma = NA))
#>   from to split_name split_value source rate_name rate_value mult_inf_probs
#> 1    S  I         NA          NA      I      <NA>         NA          FALSE
#> 2    I  R         NA          NA   NULL     gamma         NA          FALSE

# Split the infectious compartment into symptomatic (I_s) and asymptomatic
# (I_a) infections, with separate intra-household infection probabilities.
make_infection_model(
  transmit(from = "S",
           to = c("I_s", "I_a"),
           source = c("I_s", "I_a"),
           split = "phi"),
  progress(from = "I_s", to = "R", gamma_s = NA),
  progress(from = "I_a", to = "R", gamma_a = NA),
  mult_inf_probs = TRUE)
#>   from  to split_name split_value   source rate_name rate_value mult_inf_probs
#> 1    S I_s        phi          NA I_s, I_a      <NA>         NA           TRUE
#> 2    S I_a      1-phi          NA I_s, I_a      <NA>         NA           TRUE
#> 3  I_s   R       <NA>          NA     NULL   gamma_s         NA           TRUE
#> 4  I_a   R       <NA>          NA     NULL   gamma_a         NA           TRUE
```
