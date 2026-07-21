# Apply a core allocation to sampler options (rstan backend)

Configures the rstan backend for an
[`optimal_alloc()`](https://accidda.github.io/hestia/reference/optimal_alloc.md)
split: the number of chains to run in parallel is passed as `cores`, and
the per-chain thread count is exported as the `STAN_NUM_THREADS`
environment variable, which rstan reads at run time.
[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md)
owns parallelism (the parallelism arguments are rejected by
[`stan_options()`](https://accidda.github.io/hestia/reference/stan_options.md)),
so `cores` is set from the allocation. Returns the updated options.

Setting `STAN_NUM_THREADS` is a side effect;
[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md)
captures and restores the previous value so a fit does not leak its
thread count into the rest of the session. When the cmdstanr backend is
adopted, this is where the allocation is instead written as native
`parallel_chains` / `threads_per_chain` arguments, branching on the
backend.

## Usage

``` r
configure_threading(stan_opts, alloc)
```

## Arguments

- stan_opts:

  a
  [`stan_options()`](https://accidda.github.io/hestia/reference/stan_options.md)
  list.

- alloc:

  an
  [`optimal_alloc()`](https://accidda.github.io/hestia/reference/optimal_alloc.md)
  result.

## Value

`stan_opts`, with `cores` set from the allocation.
