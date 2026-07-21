# Split available cores between chain-parallelism and within-chain threads

Given the number of MCMC `chains` and the cores available, decide how
many chains to run in parallel and how many threads each chain gets for
the `reduce_sum` forward algorithm. Chain-parallelism has no
`reduce_sum` overhead and scales ~linearly, so it is filled first; any
leftover cores become per-chain threads.
[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md)
uses this to configure threading automatically for either backend.

## Usage

``` r
optimal_alloc(chains, cores = NULL)
```

## Arguments

- chains:

  number of MCMC chains (a single positive integer).

- cores:

  total cores to use. If `NULL` (default), the cores *available to the
  process* are detected with
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html),
  which respects the HPC scheduler's allocation (`SLURM_CPUS_PER_TASK`,
  PBS, SGE, LSF), cgroup CPU quotas, `getOption("mc.cores")`, and
  returns 2 under `R CMD check`. This is deliberately *not*
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html),
  which reports the whole node and would over-subscribe a scheduled HPC
  job. An explicit `cores` is used as given (no cap).

## Value

a list with integer elements `parallel_chains` and `threads_per_chain`.

## Details

`parallel_chains = min(chains, cores)` and
`threads_per_chain = floor(cores / parallel_chains)` (at least 1). So
threading only "turns on" (more than one thread) when there are more
cores than chains; otherwise every core goes to running chains in
parallel. Cores that don't divide evenly are left idle rather than
rebalanced.
