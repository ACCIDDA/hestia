test_that("stan_options() defaults chains and is otherwise empty", {
  expect_identical(stan_options(), list(chains = 4L))
  expect_identical(stan_options()$chains, 4L)
})

test_that("stan_options() passes through and coerces sampler arguments", {
  opts <- stan_options(iter = 500, chains = 2, thin = 3)

  expect_type(opts, "list")
  expect_identical(opts$iter, 500L)
  expect_identical(opts$chains, 2L)
  expect_identical(opts$thin, 3)
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

test_that("stan_options() rejects parallelism arguments (owned by run_model)", {
  expect_error(stan_options(cores = 4), "threading")
  expect_error(stan_options(parallel_chains = 4), "threading")
  expect_error(stan_options(threads_per_chain = 2), "threading")
})

test_that("stan_options() rejects non-positive or non-scalar integer arguments", {
  expect_error(stan_options(iter = 0), "iter")
  expect_error(stan_options(chains = -1), "chains")
  expect_error(stan_options(warmup = c(1, 2)), "warmup")
})

test_that("stan_options() validates seed", {
  expect_identical(stan_options(seed = 123)$seed, 123L)
  expect_error(stan_options(seed = c(1, 2)), "seed")
  expect_error(stan_options(seed = "abc"), "seed")
})
