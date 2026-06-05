# Define a infection transition

Defines a state transition in the infection process model which is the
result of a transmission (infection) event

## Usage

``` r
transmit(from, to, source = NA, split = NA)
```

## Arguments

- from:

  string giving name of origin compartment

- to:

  string giving the name of the destination compartment

- source:

  string (or vector of strings) designating which compartments are
  infectious. If NA, the destination compartment is presumed to be the
  infectious compartment.

- split:

  an optional character or numeric vector indicating what proportion of
  individuals transition into each of the `to` compartments

## Value

Data frame with a row for each unique transition and the following
columns:

- "from", the name of the origin compartment

- "to", the name of the destination compartment

- "split_name", name of split parameter if fitting

- "split_value", numeric value of split if not fitting

- "source", names of compartment(s) that are the source of new
  infections i.e. the infectious compartments

## Examples

``` r
# Specify S > E transition in an SEIR model
transmit(from = "S", to = "E", source = "I")
#>   from to split_name split_value source
#> 1    S  E         NA          NA      I

# Upon infection, individuals either enter Is (symptomatic infections) or
# Ia (asymptomatic infection) where 30% of infections are symptomatic.
transmit(from = "S", to = c("Is", "Ia"), split = 0.3)
#>   from to split_name split_value source
#> 1    S Is         NA         0.3 Is, Ia
#> 2    S Ia         NA         0.7 Is, Ia
```
