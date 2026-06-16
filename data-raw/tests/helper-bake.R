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

# Loose tolerance for the Tier-2 golden-means regression. The golden is a small
# stored vector of posterior means from the same seeded budget; the comparison is
# means only (no draws) so it stays backend-agnostic for a future rstan ->
# cmdstan swap.
golden_tol <- 0.05

# Fit a recipe with the tiny test budget. `spec` is a build_*_res_spec() list;
# `cov` toggles the covariate path (passes ih_cov/eh_cov sliced positionally).
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
    # sir_cov$covariates is positional (no hh_id): one row per individual in
    # (hh_id, part_id) order, so households 1..3 are the first
    # length(unique(part)) rows of the covariate frame.
    n_ind <- length(unique(paste(data_sub$hh_id, data_sub$part_id)))
    cov_sub <- spec$ih_cov[seq_len(n_ind), , drop = FALSE]
    args$ih_cov <- cov_sub
    args$eh_cov <- cov_sub
  }

  suppressWarnings(suppressMessages(do.call(run_model, args)))
}

# Posterior means as a named numeric vector, in a stable variable order. This is
# what both the stored golden and the live fit are reduced to for Tier 2.
draws_means <- function(draws) {
  mat <- posterior::as_draws_matrix(draws)
  means <- colMeans(mat)
  means[order(names(means))]
}

# Load the stored golden means for one recipe (NULL if not yet generated).
read_golden <- function(name) {
  path <- tryCatch(
    data_raw_path("tests", "golden-means.rds"),
    error = function(e) NULL
  )
  if (is.null(path)) {
    return(NULL)
  }
  readRDS(path)[[name]]
}
