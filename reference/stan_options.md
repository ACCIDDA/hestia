# Stan sampler options for `run_model()`

Collects and validates the arguments passed through to the Stan sampler.

Any argument accepted by
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
may be supplied, with the exception of `object` and `data` (which
[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md)
constructs internally) and `init` (which is passed to
[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md)
directly, as its defaults depend on the model structure). Sampler
controls such as `adapt_delta` and `max_treedepth` are set through
`control = list(...)`, exactly as in
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## Usage

``` r
stan_options(...)
```

## Arguments

- ...:

  arguments forwarded to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html),
  for example `iter`, `chains`, `cores`, `seed`, or
  `control = list(adapt_delta = 0.95, max_treedepth = 12)`.

## Value

a named list of validated arguments for
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## See also

[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md),
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

## Examples

``` r
stan_options()
#> $chains
#> [1] 4
#> 
stan_options(chains = 2, iter = 500)
#> $chains
#> [1] 2
#> 
#> $iter
#> [1] 500
#> 
stan_options(control = list(adapt_delta = 0.95, max_treedepth = 12))
#> $control
#> $control$adapt_delta
#> [1] 0.95
#> 
#> $control$max_treedepth
#> [1] 12
#> 
#> 
#> $chains
#> [1] 4
#> 
```
