# Define a non-infection transition

Defines a state transition in the infection process model which does not
represent an transmission (infection) event

## Usage

``` r
progress(from, to, split = NA, ...)
```

## Arguments

- from:

  string giving name of origin compartment

- to:

  string or vector of strings giving the name of the destination
  compartment

- split:

  an optional character or numeric vector indicating what proportion of
  individuals transition into each of the `to` compartments

- ...:

  Argument specifying the name and optionally the value of the
  transition rate parameter. For example if the parameter name is
  "gamma" and it is being fit this argument should be `gamma = NA`. If
  we want to provide a numeric value, x, for "gamma" this argument
  should be `gamma = x`.

## Value

Data frame with a row for each unique transition and the following
columns:

- "from", the name of the origin compartment

- "to", the name of the destination compartment

- "source", NA (so output columns match with `transmit`)

- "rate_name", name of transition rate parameter if fitting

- "rate_value", numeric value of transition rate if not fitting

- "split_name", name of split parameter if fitting

- "split_value", numeric value of split if not fitting

## Examples

``` r
# Transition from compartment A to B which occurs at rate gamma where
# gamma is a user-specified value of 0.2
progress(from = "A", to = "B", gamma = 0.2)
#>   from to rate_name rate_value split_name split_value
#> 1    A  B     gamma        0.2         NA          NA

# Transition from compartment A to B which occurs at rate gamma where
# gamma is fit
progress(from = "A", to = "B", gamma = NA)
#>   from to rate_name rate_value split_name split_value
#> 1    A  B     gamma         NA         NA          NA

# Transition from compartment A split into two destination compartments, B1
# and B2, were individuals leave A at rate delta (fit) and phi (fit) is the
# proportion of those who leave A who go into B1 (i.e. 1-phi go into B2).
# gamma is fit
progress(from = "A", to = c("B1", "B2"), split = "phi", delta = NA)
#>   from to rate_name rate_value split_name split_value
#> 1    A B1     delta         NA        phi          NA
#> 2    A B2     delta         NA      1-phi          NA

```
