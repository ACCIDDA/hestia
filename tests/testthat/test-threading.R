test_that("optimal_alloc fills chains first, then hands leftovers to threads", {
  expect_equal(
    optimal_alloc(4, cores = 8),
    list(parallel_chains = 4L, threads_per_chain = 2L)
  )
  expect_equal(
    optimal_alloc(2, cores = 8),
    list(parallel_chains = 2L, threads_per_chain = 4L)
  )
  expect_equal(
    optimal_alloc(4, cores = 16),
    list(parallel_chains = 4L, threads_per_chain = 4L)
  )
  expect_equal(
    optimal_alloc(1, cores = 8),
    list(parallel_chains = 1L, threads_per_chain = 8L)
  )
})

test_that("optimal_alloc gives one thread when there are no spare cores", {
  expect_equal(
    optimal_alloc(4, cores = 4),
    list(parallel_chains = 4L, threads_per_chain = 1L)
  )
  # more chains than cores: chains queue, still no threading
  expect_equal(
    optimal_alloc(8, cores = 4),
    list(parallel_chains = 4L, threads_per_chain = 1L)
  )
})

test_that("optimal_alloc floors (leaves idle cores) when cores don't divide", {
  expect_equal(
    optimal_alloc(4, cores = 6),
    list(parallel_chains = 4L, threads_per_chain = 1L)
  )
  expect_equal(
    optimal_alloc(3, cores = 10),
    list(parallel_chains = 3L, threads_per_chain = 3L)
  )
})

test_that("explicit cores are used as given, bypassing availableCores()", {
  # explicit cores skip parallelly entirely -- deterministic regardless of the
  # machine, the scheduler, or R CMD check.
  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_equal(
    optimal_alloc(4, cores = 8),
    list(parallel_chains = 4L, threads_per_chain = 2L)
  )
})

test_that("auto-detected cores respect availableCores() under R CMD check", {
  # parallelly::availableCores() returns 2 when _R_CHECK_LIMIT_CORES_ is set.
  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_equal(
    optimal_alloc(4),
    list(parallel_chains = 2L, threads_per_chain = 1L)
  )
})

test_that("auto-detected cores respect a constraining mc.cores", {
  # availableCores() takes the *minimum* across sources, so mc.cores = 1 is a
  # hard floor independent of the test machine's core count.
  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "")
  withr::local_options(mc.cores = 1)
  expect_equal(
    optimal_alloc(4),
    list(parallel_chains = 1L, threads_per_chain = 1L)
  )
})

test_that("the auto-detect path returns a valid allocation", {
  alloc <- optimal_alloc(4)
  expect_named(alloc, c("parallel_chains", "threads_per_chain"))
  expect_true(alloc$parallel_chains >= 1L && alloc$parallel_chains <= 4L)
  expect_true(alloc$threads_per_chain >= 1L)
})

test_that("optimal_alloc validates chains", {
  expect_error(optimal_alloc(0, cores = 8))
  expect_error(optimal_alloc(-1, cores = 8))
  expect_error(optimal_alloc(2.5, cores = 8))
  expect_error(optimal_alloc(c(2, 4), cores = 8))
})

test_that("configure_threading sets cores and STAN_NUM_THREADS from the alloc", {
  withr::local_envvar(STAN_NUM_THREADS = "")
  opts <- configure_threading(
    list(chains = 4L),
    list(parallel_chains = 4L, threads_per_chain = 2L)
  )
  expect_equal(opts$cores, 4L)
  expect_equal(Sys.getenv("STAN_NUM_THREADS"), "2")
})

test_that("configure_threading owns cores, overwriting any stray value", {
  # run_model owns parallelism (stan_options() rejects `cores`), so the
  # allocation's parallel_chains wins even if a raw list smuggles one in.
  withr::local_envvar(STAN_NUM_THREADS = "")
  opts <- configure_threading(
    list(chains = 4L, cores = 8L),
    list(parallel_chains = 4L, threads_per_chain = 2L)
  )
  expect_equal(opts$cores, 4L)
  expect_equal(Sys.getenv("STAN_NUM_THREADS"), "2")
})

test_that("configure_threading exports a single thread when threading is off", {
  withr::local_envvar(STAN_NUM_THREADS = "")
  opts <- configure_threading(
    list(chains = 4L),
    list(parallel_chains = 4L, threads_per_chain = 1L)
  )
  expect_equal(Sys.getenv("STAN_NUM_THREADS"), "1")
})
