# Backend handling for run_model(): the cross-backend vocabulary guard and the
# per-backend fit functions. Each fit function forwards the user's stan_options()
# verbatim to that backend's native sampler, so calls feel like using rstan /
# cmdstanr directly.

# cmdstanr-only argument names, shown when the active backend is rstan (i.e. the
# user reached for a cmdstanr word), each mapped to the rstan way to do it.
cmdstanr_hints <- c(
  parallel_chains   = "use `cores`",
  iter_warmup       = "use `iter` (with `warmup`)",
  iter_sampling     = "use `iter` (with `warmup`)",
  adapt_delta       = "set inside `control = list(adapt_delta = ...)`",
  max_treedepth     = "set inside `control = list(max_treedepth = ...)`",
  step_size         = "set inside `control = list(stepsize = ...)`",
  threads_per_chain = "rstan parallelises chains via `cores`",
  output_dir        = "no rstan equivalent",
  sig_figs          = "no rstan equivalent"
)

# rstan-only argument names, shown when the active backend is cmdstanr, mapped to
# the cmdstanr way to do it.
rstan_hints <- c(
  cores           = "use `parallel_chains`",
  control         = "set `adapt_delta`/`max_treedepth`/`step_size` as top-level arguments",
  iter            = "use `iter_warmup` and `iter_sampling`",
  warmup          = "use `iter_warmup`",
  pars            = "not supported; use `save_states` in run_model()",
  include         = "not supported; use `save_states` in run_model()",
  sample_file     = "no cmdstanr equivalent",
  diagnostic_file = "no cmdstanr equivalent"
)

#' Reject the other backend's argument vocabulary
#'
#' @param arg_names names of the arguments supplied to [stan_options()].
#' @param backend the backend the options are being built for.
#' @keywords internal
check_backend_vocab <- function(arg_names, backend) {
  foreign <- switch(
    backend,
    rstan    = cmdstanr_hints,
    cmdstanr = rstan_hints
  )
  bad <- intersect(arg_names, names(foreign))
  if (length(bad) > 0) {
    bullets <- paste0("  - `", bad, "`: ", foreign[bad], collapse = "\n")
    stop(
      "These stan_options() arguments are not valid for the '", backend,
      "' backend:\n", bullets,
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Error if the selected backend's package is not installed
#'
#' rstan is always available (a hard dependency); cmdstanr is optional, so
#' selecting it without the package installed fails early here rather than deep
#' inside the fit.
#'
#' @param backend the backend to check.
#' @keywords internal
check_backend_available <- function(backend) {
  if (backend == "cmdstanr" && !requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "backend = 'cmdstanr' requires the cmdstanr package, which is not ",
      "installed. Install it from https://mc-stan.org/cmdstanr/, or use ",
      "backend = 'rstan'.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Dispatch a fit to the chosen backend
#'
#' @param backend one of `"rstan"` or `"cmdstanr"`.
#' @param model_name `"hmm"` or `"hmm_cov"`.
#' @param dat_stan the Stan data list from `make_stan_data()`.
#' @param init the init list, sized to the chain count.
#' @param stan_opts the validated, backend-tagged `stan_options()` list.
#' @param save_states whether to keep the per-timestep state probabilities.
#' @returns the backend's fit object (a `stanfit` or `CmdStanMCMC`).
#' @keywords internal
fit_model <- function(backend, model_name, dat_stan, init, stan_opts, save_states) {
  switch(
    backend,
    rstan    = fit_rstan(model_name, dat_stan, init, stan_opts, save_states),
    cmdstanr = fit_cmdstanr(model_name, dat_stan, init, stan_opts, save_states),
    stop("Unknown backend: ", backend, call. = FALSE)
  )
}

#' @keywords internal
fit_rstan <- function(model_name, dat_stan, init, stan_opts, save_states) {
  args <- stan_opts
  attr(args, "hestia_backend") <- NULL
  args$object <- stanmodels[[model_name]]
  args$data   <- dat_stan
  args$init   <- init
  if (!save_states) {
    # Drop the per-timestep state probabilities from the saved output.
    args$pars    <- "logalpha"
    args$include <- FALSE
  }
  do.call(rstan::sampling, args)
}

#' @keywords internal
fit_cmdstanr <- function(model_name, dat_stan, init, stan_opts, save_states) {
  # nocov start: needs the CmdStan toolchain, unavailable on CI/CRAN
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "The 'cmdstanr' backend requires the cmdstanr package. ",
      "See https://mc-stan.org/cmdstanr/#installing-the-r-package, ",
      "or use backend = 'rstan'.",
      call. = FALSE
    )
  }
  if (!save_states) {
    warning(
      "save_states is not supported by the cmdstanr backend; ",
      "logalpha is always written to the output.",
      call. = FALSE
    )
  }
  stan_file <- system.file(
    "stan", paste0(model_name, ".stan"),
    package = "hestia", mustWork = TRUE
  )
  # Compile into a writable user cache, not next to the installed .stan file
  # (the package directory may be read-only, and stray executables there trip
  # R CMD check's "executable files" warning). cmdstan_model() reuses the cached
  # executable across sessions and recompiles only when the .stan source is
  # newer than it -- e.g. after a package update reinstalls the .stan file.
  exe_dir <- tools::R_user_dir("hestia", "cache")
  dir.create(exe_dir, showWarnings = FALSE, recursive = TRUE)
  mod <- cmdstanr::cmdstan_model(stan_file, dir = exe_dir)

  args <- stan_opts
  attr(args, "hestia_backend") <- NULL
  args$data <- dat_stan
  args$init <- init
  do.call(mod$sample, args)
  # nocov end
}
