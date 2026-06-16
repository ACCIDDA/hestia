# Shared setup for the data-raw bake integration tests.
#
# These tests run ONE tiny seeded Stan fit per bundled-dataset recipe and guard
# the data-raw/*_res.R bake scripts against silent API breakage. They live under
# data-raw/ (build-excluded) and run via the dedicated data-raw-integration
# workflow, not R CMD check, so they never slow the package check or covr.
#
# The model definitions come from data-raw/recipes.R, the SAME file the
# production bakes source, so a #57-class change to run_model()'s output fails
# both the bake and these tests. Only the sampler budget differs: the bakes use
# stan_options(iter = 1000, cores = 4); the tests use test_stan_opts below (one
# short, seeded chain on a ~3-household slice).
#
# pkgload::load_all() puts the package's exported and internal objects (run_model,
# the make_* constructors, the bundled data) on the search path, so recipes.R and
# these helpers can call them directly without hestia::.

# Resolve a path that lives under data-raw/, robust to the caller's working
# directory. testthat::test_dir() runs helpers with the working directory set to
# the test folder (data-raw/tests/), while the bake and regen scripts run from
# the package root (data-raw/...); try both layouts.
data_raw_path <- function(...) {
  rel <- file.path(...)
  candidates <- c(
    file.path("data-raw", rel),
    file.path("..", rel),
    rel
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("could not locate data-raw file: ", rel, call. = FALSE)
  }
  hit[[1L]]
}

source(data_raw_path("recipes.R"))

# Tiny, fully pinned sampler budget. chains = 1 keeps the fit cheap and, with a
# fixed seed, deterministic: run_model() sizes its init list to stan_opts$chains
# (R/hestia_functions.R), so a single-chain fit goes straight through run_model()
# the same way the production bake does.
test_stan_opts <- stan_options(chains = 1, iter = 80, seed = 1L, refresh = 0)

# Keep the slice to households 1..3. make_stan_data() uses max(hh_id) as the
# household count and groups by hh_id, so the slice must preserve contiguous
# 1..N ids; hh_id <= 3 does.
test_n_hh <- 3L

# Per-parameter RELATIVE tolerance (with a small absolute floor) for the Tier-2
# golden-means regression. waldo's pooled tolerance is scaled by the whole
# vector's magnitude, so a large exp-scale coefficient would swamp the bound on a
# small probability; expect_means_close() below compares each parameter on its
# own so the same headroom applies to every variable.
golden_tol <- 0.05
golden_floor <- 1e-3

# The two heavier fits (siir, sir_cov) only run when HESTIA_RUN_INTEGRATION is
# set, so a bare local test_dir() pays for just the cheap SIR fit; the dedicated
# CI workflow sets it so all three run there. data-raw is build-excluded, so
# these tests never reach R CMD check / CRAN: skip_on_cran() would guard a
# context that cannot occur.
run_full_integration <- function() {
  isTRUE(as.logical(Sys.getenv("HESTIA_RUN_INTEGRATION", "false")))
}

# Slice a positional covariate frame to the first n_ind individuals. The frame
# carries one row per individual in the full data, in (hh_id, part_id) order;
# assert that one-row-per-individual shape so a reordered or regenerated dataset
# fails loudly here rather than silently mis-aligning the slice.
slice_cov <- function(cov, data, n_ind) {
  n_full <- length(unique(paste(data$hh_id, data$part_id)))
  stopifnot(nrow(cov) == n_full)
  cov[seq_len(n_ind), , drop = FALSE]
}

# Fit a recipe with the tiny test budget. `spec` is a build_*_res_spec() list;
# `cov` toggles the covariate path (slices and passes ih_cov and eh_cov
# independently from the spec, so a recipe that diverged them would be caught).
fit_test_recipe <- function(spec, cov = FALSE) {
  data_sub <- spec$data[spec$data$hh_id <= test_n_hh, ]

  args <- list(
    inf_model = spec$inf_model,
    obs_model = spec$obs_model,
    data = data_sub,
    init_probs = spec$init_probs,
    stan_opts = test_stan_opts
  )

  if (cov) {
    n_ind <- length(unique(paste(data_sub$hh_id, data_sub$part_id)))
    args$ih_cov <- slice_cov(spec$ih_cov, spec$data, n_ind)
    args$eh_cov <- slice_cov(spec$eh_cov, spec$data, n_ind)
  }

  suppressWarnings(suppressMessages(do.call(run_model, args)))
}

# Posterior means as a named numeric vector in a stable variable order. This is
# what both the stored golden and the live fit reduce to for Tier 2.
draws_means <- function(draws) {
  summ <- posterior::summarise_draws(draws, mean)
  out <- stats::setNames(summ$mean, summ$variable)
  out[order(names(out))]
}

# Load the stored golden means for one recipe.
read_golden <- function(name) {
  readRDS(data_raw_path("tests", "golden-means.rds"))[[name]]
}

# Tier-2 regression assertion: every parameter's posterior mean is within a
# relative tolerance (with an absolute floor) of the stored golden. Per-parameter
# so the same headroom applies to a small probability and a large coefficient.
expect_means_close <- function(means, golden) {
  testthat::expect_setequal(names(means), names(golden))
  golden <- golden[names(means)]
  bound <- pmax(abs(golden) * golden_tol, golden_floor)
  off <- abs(means - golden) > bound
  testthat::expect_false(
    any(off),
    info = paste(
      "means outside tolerance:",
      paste(names(means)[off], collapse = ", ")
    )
  )
}
