#' @title Stan sampler options for `run_model()`
#'
#' @description
#' Collects and validates sampler arguments for the chosen `backend`, forwarding
#' them **verbatim** so calls feel native to that backend. Use the backend's own
#' argument names; mixing one backend's vocabulary into the other errors with a
#' hint. `object`, `data`, and `init` may not be set here — [run_model()] owns
#' those.
#'
#' @param ... sampler arguments forwarded verbatim to the chosen backend's
#'   sampler. Use the backend's own names: for `"rstan"`, the
#'   [rstan::sampling()] arguments (`iter`, `cores`, `seed`,
#'   `control = list(adapt_delta = 0.95)`); for `"cmdstanr"`, the `$sample()`
#'   arguments (`iter_warmup`, `iter_sampling`, `parallel_chains`,
#'   `adapt_delta`, ...).
#' @param backend which Stan interface to target, one of `"rstan"` (default) or
#'   `"cmdstanr"`. Determines which argument vocabulary is accepted and which
#'   sampler [run_model()] calls. Whatever is set here wins. Selecting
#'   `"cmdstanr"` errors if the cmdstanr package is not installed.
#' @param chains number of MCMC chains. Native to both backends; also used by
#'   [run_model()] to size the per-chain list of initial values.
#'
#' @returns a named list of validated sampler arguments, tagged with the
#'   backend it was built for.
#'
#' @examples
#' stan_options()
#' stan_options(chains = 2, iter = 500)
#' stan_options(control = list(adapt_delta = 0.95, max_treedepth = 12))
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   stan_options(backend = "cmdstanr", parallel_chains = 4, iter_warmup = 500)
#' }
#'
#' @seealso [run_model()], [rstan::sampling()]
#' @export
stan_options <- function(..., backend = "rstan", chains = 4L) {
  backend <- check_backend_available(backend)
  res <- list(...)

  if ("object" %in% names(res)) {
    stop(
      "Passing 'object' in stan_options() is not allowed; ",
      "The model object should be passed in `run_model()` instead.",
      call. = FALSE
    )
  }
  if ("data" %in% names(res)) {
    stop(
      "Passing 'data' in stan_options() is not allowed; ",
      "run_model() constructs the 'data' argument internally.",
      call. = FALSE
    )
  }
  if ("init" %in% names(res)) {
    stop(
      "Passing 'init' in stan_options() is not allowed; ",
      "pass it to run_model() directly, as its defaults depend on the ",
      "model structure.",
      call. = FALSE
    )
  }

  # Reject the other backend's vocabulary with a "did you mean" hint.
  check_backend_vocab(names(res), backend)

  # `chains` is native to both backends and run_model() needs it to size the
  # init list, so it is an explicit argument with a default rather than a
  # value injected from `...`.
  res$chains <- check_positive_int(chains, "chains")

  # Validate the remaining positive-integer count arguments native to this
  # backend.
  for (arg in intersect(names(res), backend_int_args(backend))) {
    res[[arg]] <- check_positive_int(res[[arg]], arg)
  }

  if ("seed" %in% names(res)) {
    res[["seed"]] <- check_seed(res[["seed"]])
  }

  attr(res, "hestia_backend") <- backend
  res
}
