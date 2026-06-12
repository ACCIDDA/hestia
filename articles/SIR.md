# Basic SIR Model

## SIR model

### Simulate data

We will use data simulated from a basic SIR model for 500 households. We
will assume that all individuals start susceptible. The per day
intra-household infection probability is 5% and the daily
extra-household infection probability is 1%.

``` r

head(sir)
```

    ##   t part_id enroll state hh_size hh_id pcr igg
    ## 1 1       1      1     1       3     1   0   1
    ## 2 1       2      1     1       3     1   0   0
    ## 3 1       3      0     1       3     1   0   0
    ## 4 2       1      1     1       3     1   0   0
    ## 5 2       2      1     1       3     1   0   0
    ## 6 2       3      0     1       3     1   0   0

### Model specification

In addition to the data, the user needs to specify the underlying model
structure, which involves two components:

1.  **The infection process model** which specifies the underlying
    compartmental model structure.

2.  **The observation process model** which specifies the probability of
    observing a positive outcome given an individual is in each
    compartment of the infection process model.

#### Infection process model

The infection process model is built up from a series of user-defined
transitions. There are two types of transitions:
[`transmit()`](https://accidda.github.io/hestia/reference/transmit.md)
defines a transition that occurs as a result of infection (e.g. S\>I)
and
[`progress()`](https://accidda.github.io/hestia/reference/progress.md)
which defines a non-infection transition that occurs at a constant rate
(e.g. I-\>R). For both functions, the user specifies a destination (e.g
`from = S`) and source (e.g `to = I`) compartment. Additionally for
[`progress()`](https://accidda.github.io/hestia/reference/progress.md)
transitions, the user must enter a parameter for the transition rate.
The name of this parameter is at the discretion of the user
(e.g. `gamma` for recovery rate) and can be set to a specific numeric
value of to NA if model should fit that parameter.

``` r

# Basic SIR, fit recovery rate
inf_process <- make_infection_model(transmit(from = "S", to = "I"),
                                    progress(from = "I", to = "R", gamma = NA))
```

#### Observation process model

The observation model is composed of a series of named vectors. Each
vector corresponds to an observation type. Each entry in the vector is
named for the corresponding compartment in the infection process model.
The value is the probability of observing a positive observation given
the individual is in the compartment.

In this example we have two observation types, one which is likely to be
positive when a person is actively infectious (e.g. PCR test results)
and the other which is more likely to be positive when the person is
recovered and immune to infection (e.g. IgG antibody).

``` r

obs_process <- make_observation_model(
  pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
  igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
)
```

Now we are ready to run the model!

``` r

# Note this is very computationally intensive so the block is not evaluated
# when knitting

## Inputs
# 1. Infection process model
# 2. Observation process model
# 3. Raw data
# 4. Covariates (optional)
# 5. Starting state probabilities

sir_res <- run_model(
  inf_model = inf_process,
  obs_model = obs_process,
  data = sir,
  init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10),
  stan_opts = stan_options(iter = 1000, cores = 4)
)
```

Now let’s look at the model results using the `posterior` package.

``` r

# The output from the model in the previous code chunk is available in the
# sir_res package data object

# Get summary statistics. Parameters are returned on the natural (model) scale,
# so no back-transformation is needed.
draws_sum <- summarise_draws(sir_res, "mean", "median", ~ quantile2(.x, probs = c(0.025, 0.975))) |>
  mutate(var_type = c(rep("Infection probability", 2), "Recovery rate"))

# "True" values from underlying simulation
draws_sum$true_value <- c(0.01, 0.05, 0.2)

draws_sum
```

    ## # A tibble: 3 × 7
    ##   variable  mean median  q2.5 q97.5 var_type              true_value
    ##   <chr>    <dbl>  <dbl> <dbl> <dbl> <chr>                      <dbl>
    ## 1 eh_prob  -4.61  -4.61 -4.68 -4.54 Infection probability       0.01
    ## 2 ih_prob  -2.85  -2.85 -2.96 -2.74 Infection probability       0.05
    ## 3 gamma    -1.37  -1.37 -1.44 -1.31 Recovery rate               0.2

### Add in covariates

We can add covariates on the intra- and extra-household infection risk
by including the optional `ih_cov` and `eh_cov` arguments, respectively,
in
[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md).

``` r

# Simulated data from a basic SIR model for 500 households
head(sir_cov$observations)

# Run model with x1 and x2 as covariates
sir_cov_res <- run_model(
  inf_model = inf_process,
  obs_model = obs_process,
  data = sir_cov$observations,
  ih_cov = sir_cov$covariates,
  eh_cov = sir_cov$covariates,
  init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10),
  stan_opts = stan_options(iter = 1000, cores = 4)
)
```

Let’s take a look at the results once more.

``` r

# The output from the model in the previous code chunk is available in the
# sir_cov_res package data object

# Get summary statistics. Infection probabilities and rates are returned on the
# natural (model) scale and coefficients on the natural (exponentiated) scale,
# so no back-transformation is needed.
draws_sum <- summarise_draws(sir_cov_res, "mean", "median",
                             ~ quantile2(.x, probs = c(0.025, 0.975))) |>
  mutate(var_type = c(rep("Infection probability", 2), "Recovery rate", rep("Coefficient", 4)))

# "True" values from underlying simulation
draws_sum$true_vals <- c(0.01, 0.05, 0.2, exp(-0.4), exp(0.7), exp(0.8), exp(0.1))

draws_sum
```

    ## # A tibble: 7 × 7
    ##   variable   mean median   q2.5   q97.5 var_type              true_vals
    ##   <chr>     <dbl>  <dbl>  <dbl>   <dbl> <chr>                     <dbl>
    ## 1 eh_prob  -4.64  -4.64  -4.75  -4.52   Infection probability     0.01 
    ## 2 ih_prob  -2.92  -2.92  -3.09  -2.75   Infection probability     0.05 
    ## 3 gamma    -1.40  -1.41  -1.47  -1.34   Recovery rate             0.2  
    ## 4 x1_eh    -0.460 -0.458 -0.620 -0.295  Coefficient               0.670
    ## 5 x2_eh     0.709  0.708  0.562  0.859  Coefficient               2.01 
    ## 6 x1_ih     0.778  0.782  0.558  0.990  Coefficient               2.23 
    ## 7 x2_ih    -0.170 -0.169 -0.420  0.0772 Coefficient               1.11
