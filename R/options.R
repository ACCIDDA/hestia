#' @title Stan sampler options for `run_model()`
#'
#' @description
#' Collects and validates the arguments passed through to the Stan sampler.
#'
#' Any argument accepted by [rstan::sampling()] may be supplied, with the
#' exception of the arguments [run_model()] manages itself: `object` and `data`
#' (constructed internally), `init` (passed to [run_model()] directly, as its
#' defaults depend on the model structure), and the parallelism controls `cores`,
#' `parallel_chains`, and `threads_per_chain` (owned by [run_model()]'s
#' `threading` argument, which splits the available cores between chains and
#' threads). Sampler controls such as `adapt_delta` and `max_treedepth` are set
#' through `control = list(...)`, exactly as in [rstan::sampling()].
#'
#' @param ... arguments forwarded to [rstan::sampling()], for example `iter`,
#'   `chains`, `seed`, or
#'   `control = list(adapt_delta = 0.95, max_treedepth = 12)`.
#'
#' @returns a named list of validated arguments for [rstan::sampling()].
#'
#' @examples
#' stan_options()
#' stan_options(chains = 2, iter = 500)
#' stan_options(control = list(adapt_delta = 0.95, max_treedepth = 12))
#'
#' @seealso [run_model()], [rstan::sampling()]
#' @export
stan_options <- function(...) {
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

  # Parallelism is owned by run_model(): it splits the available cores between
  # chains and threads for the active backend (see its `threading` argument), so
  # these must not be set here.
  parallel_args <- c("cores", "parallel_chains", "threads_per_chain")
  bad <- intersect(names(res), parallel_args)
  if (length(bad) > 0) {
    stop(
      "Passing ", paste0("'", bad, "'", collapse = ", "),
      " in stan_options() is not allowed; parallelism is managed by ",
      "run_model() -- use its `threading` argument.",
      call. = FALSE
    )
  }

  if (is.null(res$chains)) {
    res$chains <- 4L
  }

  int_args <- c("iter", "chains", "warmup")
  for (arg in intersect(names(res), int_args)) {
    res[[arg]] <- check_positive_int(res[[arg]], arg)
  }

  if ("seed" %in% names(res)) {
    res[["seed"]] <- check_seed(res[["seed"]])
  }

  res
}
