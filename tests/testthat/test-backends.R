# Backend selection + cross-backend vocabulary guard. These exercise
# stan_options() and check_backend_vocab() only -- no Stan is fit.

test_that("cmdstanr backend accepts its own native arguments", {
  opts <- stan_options(
    backend = "cmdstanr",
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.95
  )

  expect_identical(attr(opts, "hestia_backend"), "cmdstanr")
  expect_identical(opts$parallel_chains, 4L)
  expect_identical(opts$iter_warmup, 500L)
})

test_that("rstan backend rejects cmdstanr vocabulary with a hint", {
  expect_error(
    stan_options(backend = "rstan", parallel_chains = 4),
    "parallel_chains"
  )
  expect_error(
    stan_options(backend = "rstan", iter_warmup = 500),
    "iter_warmup"
  )
})

test_that("cmdstanr backend rejects rstan vocabulary with a hint", {
  expect_error(
    stan_options(backend = "cmdstanr", cores = 4),
    "cores"
  )
  expect_error(
    stan_options(backend = "cmdstanr", control = list(adapt_delta = 0.95)),
    "control"
  )
  expect_error(
    stan_options(backend = "cmdstanr", iter = 1000),
    "iter"
  )
})

test_that("check_backend_vocab passes clean argument sets", {
  expect_silent(check_backend_vocab(c("iter", "chains", "cores"), "rstan"))
  expect_silent(
    check_backend_vocab(c("iter_warmup", "parallel_chains"), "cmdstanr")
  )
})

test_that("fit_model errors on an unknown backend", {
  expect_error(
    fit_model("nonsense", "hmm", list(), NULL, stan_options(), FALSE),
    "Unknown backend"
  )
})
