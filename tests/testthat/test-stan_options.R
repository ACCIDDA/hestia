test_that("stan_options() defaults chains and is otherwise empty", {
  opts <- stan_options()
  expect_identical(opts$chains, 4L)
  expect_named(opts, "chains")
})

test_that("stan_options() defaults to and tags the rstan backend", {
  expect_identical(attr(stan_options(), "hestia_backend"), "rstan")
  expect_identical(
    attr(stan_options(backend = "cmdstanr"), "hestia_backend"),
    "cmdstanr"
  )
  expect_error(stan_options(backend = "nonsense"))
})

test_that("stan_options() passes through and coerces sampler arguments", {
  opts <- stan_options(iter = 500, chains = 2, cores = 4)

  expect_type(opts, "list")
  expect_identical(opts$iter, 500L)
  expect_identical(opts$chains, 2L)
  expect_identical(opts$cores, 4L)
})

test_that("stan_options() forwards arbitrary sampling() arguments untouched", {
  opts <- stan_options(control = list(adapt_delta = 0.95, max_treedepth = 12))

  expect_identical(opts$control, list(adapt_delta = 0.95, max_treedepth = 12))
})

test_that("stan_options() rejects internally-managed arguments", {
  expect_error(stan_options(object = "hmm"), "object")
  expect_error(stan_options(data = list()), "data")
  expect_error(stan_options(init = list()), "init")
})

test_that("stan_options() rejects non-positive or non-scalar integer arguments", {
  expect_error(stan_options(iter = 0), "iter")
  expect_error(stan_options(chains = -1), "chains")
  expect_error(stan_options(cores = c(1, 2)), "cores")
})

test_that("stan_options() validates seed", {
  expect_identical(stan_options(seed = 123)$seed, 123L)
  expect_error(stan_options(seed = c(1, 2)), "seed")
  expect_error(stan_options(seed = "abc"), "seed")
})
