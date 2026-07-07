#' Split available cores between chain-parallelism and within-chain threads
#'
#' @description
#' Given the number of MCMC `chains` and the cores available, decide how many
#' chains to run in parallel and how many threads each chain gets for the
#' `reduce_sum` forward algorithm. Chain-parallelism has no `reduce_sum`
#' overhead and scales ~linearly, so it is filled first; any leftover cores
#' become per-chain threads. `run_model()` uses this to configure threading
#' automatically for either backend.
#'
#' @details
#' `parallel_chains = min(chains, cores)` and
#' `threads_per_chain = floor(cores / parallel_chains)` (at least 1). So
#' threading only "turns on" (more than one thread) when there are more cores
#' than chains; otherwise every core goes to running chains in parallel. Cores
#' that don't divide evenly are left idle rather than rebalanced.
#'
#' @param chains number of MCMC chains (a single positive integer).
#' @param cores total cores to use. If `NULL` (default), the cores *available to
#'   the process* are detected with [parallelly::availableCores()], which
#'   respects the HPC scheduler's allocation (`SLURM_CPUS_PER_TASK`, PBS, SGE,
#'   LSF), cgroup CPU quotas, `getOption("mc.cores")`, and returns 2 under
#'   `R CMD check`. This is deliberately *not* `parallel::detectCores()`, which
#'   reports the whole node and would over-subscribe a scheduled HPC job. An
#'   explicit `cores` is used as given (no cap).
#' @returns a list with integer elements `parallel_chains` and
#'   `threads_per_chain`.
#' @keywords internal
optimal_alloc <- function(chains, cores = NULL) {
  checkmate::assert_count(chains, positive = TRUE)
  chains <- as.integer(chains)

  if (is.null(cores)) {
    # Cores the process is *allowed* to use (scheduler-, cgroup-, and
    # check-aware), not the node's physical core count.
    cores <- parallelly::availableCores()
  }
  checkmate::assert_count(cores, positive = TRUE)
  cores <- as.integer(cores)

  parallel_chains   <- min(chains, cores)
  threads_per_chain <- max(1L, cores %/% parallel_chains)

  list(parallel_chains = parallel_chains, threads_per_chain = threads_per_chain)
}

#' Apply a core allocation to sampler options (rstan backend)
#'
#' @description
#' Configures the rstan backend for an [optimal_alloc()] split: the number of
#' chains to run in parallel is passed as `cores`, and the per-chain thread count
#' is exported as the `STAN_NUM_THREADS` environment variable, which rstan reads
#' at run time. `run_model()` owns parallelism (the parallelism arguments are
#' rejected by [stan_options()]), so `cores` is set from the allocation. Returns
#' the updated options.
#'
#' Setting `STAN_NUM_THREADS` is a side effect; `run_model()` captures and
#' restores the previous value so a fit does not leak its thread count into the
#' rest of the session. When the cmdstanr backend is adopted, this is where the
#' allocation is instead written as native `parallel_chains` / `threads_per_chain`
#' arguments, branching on the backend.
#'
#' @param stan_opts a [stan_options()] list.
#' @param alloc an [optimal_alloc()] result.
#' @returns `stan_opts`, with `cores` set from the allocation.
#' @keywords internal
configure_threading <- function(stan_opts, alloc) {
  stan_opts$cores <- alloc$parallel_chains
  Sys.setenv(STAN_NUM_THREADS = alloc$threads_per_chain)
  stan_opts
}
