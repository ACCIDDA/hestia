# Backend selection + cross-backend vocabulary guard. The vocabulary logic is
# tested directly via check_backend_vocab() so it runs without cmdstanr
# installed; the stan_options() integration tests for cmdstanr are skipped when
# cmdstanr is unavailable (e.g. on CI / CRAN).

test_that("check_backend_vocab rejects the other backend's vocabulary", {
  expect_error(check_backend_vocab("parallel_chains", "rstan"), "parallel_chains")
  expect_error(check_backend_vocab("iter_warmup", "rstan"), "iter_warmup")
  expect_error(check_backend_vocab("cores", "cmdstanr"), "cores")
  expect_error(check_backend_vocab("control", "cmdstanr"), "control")
  expect_error(check_backend_vocab("iter", "cmdstanr"), "iter")
})

test_that("check_backend_vocab passes clean argument sets", {
  expect_silent(check_backend_vocab(c("iter", "cores"), "rstan"))
  expect_silent(
    check_backend_vocab(c("iter_warmup", "parallel_chains"), "cmdstanr")
  )
})

test_that("rstan backend rejects cmdstanr vocabulary via stan_options", {
  expect_error(stan_options(backend = "rstan", parallel_chains = 4), "parallel_chains")
})

test_that("stan_options errors for cmdstanr when it is not installed", {
  skip_if(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr is installed")
  expect_error(stan_options(backend = "cmdstanr"), "cmdstanr")
})

test_that("cmdstanr backend builds and tags native options when available", {
  skip_if_not_installed("cmdstanr")
  opts <- stan_options(
    backend = "cmdstanr",
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500
  )
  expect_identical(attr(opts, "hestia_backend"), "cmdstanr")
  expect_identical(opts$parallel_chains, 4L)
  expect_identical(opts$iter_warmup, 500L)
})

test_that("cmdstanr backend rejects rstan vocabulary when available", {
  skip_if_not_installed("cmdstanr")
  expect_error(stan_options(backend = "cmdstanr", cores = 4), "cores")
  expect_error(stan_options(backend = "cmdstanr", iter = 1000), "iter")
})

test_that("fit_model errors on an unknown backend", {
  expect_error(
    fit_model("nonsense", "hmm", list(), NULL, stan_options(), FALSE),
    "Unknown backend"
  )
})
