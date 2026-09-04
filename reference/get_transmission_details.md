# Create Transmission Probability Matrix

Create Transmission Probability Matrix

## Usage

``` r
get_transmission_details(inf_model)
```

## Arguments

- inf_model:

  infection process model object yielded by make_infection_model()

## Value

A list with the following entries:

- "states", a character vector of state names

- "trans_matrix", a matrix with transition rate values between origin
  (column) and destination (row) compartments where given. Where no
  values are provided for transition rates the entry is equal to 1e-10.

- "mult_matrix", a matrix with split values for the transition between
  origin (column) and destination (row) compartments where given. Where
  no split values are provided for transition rates the entry is equal
  to 1.

- "trans_to_fit", a data frame with information on all transitions
  within the model, with parameter indexing details if the transition
  rate is being fit

- "mult_to_fit", a data frame with information on all splits specified
  in the model, with parameter indexing details if the split is being
  fit

- "inf_states", numeric values corresponding to the indices in "states"
  which are the infectious states

- "mult_ih_inf_probs", indicator for whether infectious states have
  shared or separate inta-household infection probabilities. If FALSE
  then all infection probabilities are shared across infectious
  compartments.

## Examples

``` r
# Split the infectious compartment into symptomatic (I_s) and asymptomatic
# (I_a) infections, with separate intra-household infection probabilities.
inf_model <- make_infection_model(
  transmit(from = "S",
           to = c("I_s", "I_a"),
           source = c("I_s", "I_a"),
           split = "phi"),
  progress(from = "I_s", to = "R", gamma_s = NA),
  progress(from = "I_a", to = "R", gamma_a = NA),
  mult_ih_inf_probs = TRUE)

get_transmission_details(inf_model)
#> $states
#> [1] "S"   "I_s" "I_a" "R"  
#> 
#> $trans_matrix
#>        from_S from_I_s from_I_a from_R
#> to_S    1e-10    1e-10    1e-10  1e-10
#> to_I_s  1e-10    1e-10    1e-10  1e-10
#> to_I_a  1e-10    1e-10    1e-10  1e-10
#> to_R    1e-10    1e-10    1e-10  1e-10
#> 
#> $mult_matrix
#>        from_S from_I_s from_I_a from_R
#> to_S        1        1        1      1
#> to_I_s      1        1        1      1
#> to_I_a      1        1        1      1
#> to_R        1        1        1      1
#> 
#> $trans_to_fit
#>   from  to trans_row trans_col rate_name param source
#> 1    S I_s         2         1      <NA>     0   2, 3
#> 2    S I_a         3         1      <NA>     0   2, 3
#> 3  I_s   R         4         2   gamma_s     1      0
#> 4  I_a   R         4         3   gamma_a     2      0
#> 
#> $mult_to_fit
#>   from  to mult_row mult_col mult_name param
#> 1    S I_s        2        1       phi     1
#> 2    S I_a        3        1     1-phi    -1
#> 
#> $inf_states
#> [1] 2 3
#> 
#> $mult_ih_inf_probs
#> [1] TRUE
#> 
#> $mult_eh_inf_probs
#> [1] FALSE
#> 
#> $compete
#> [1] 0 0 0 0
#> 

```
