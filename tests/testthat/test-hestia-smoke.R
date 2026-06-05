# End-to-end smoke test. This fits a real (if tiny) Stan model on the bundled
# `sir` data and checks that the pipeline returns sensible, well-named draws.
# It deliberately uses a single short chain so it is cheap enough to run in CI
# on every push (no skip_on_cran()).

test_that("the package and its bundled data are available", {
  expect_true(requireNamespace("hestia", quietly = TRUE))

  data(sir, package = "hestia", envir = environment())
  expect_s3_class(sir, "data.frame")
  expect_true(all(c("t", "part_id", "hh_size", "hh_id", "pcr", "igg") %in%
    names(sir)))
})

test_that("run_model fits an SIR model end-to-end and returns named draws", {
  inf_mod <- sir_infection_model()
  obs_mod <- sir_observation_model()
  dat <- sir_subset(10)

  suppressWarnings(suppressMessages(
    draws <- run_model(
      inf_model = inf_mod,
      obs_model = obs_mod,
      data = dat,
      init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10),
      iter = 100,
      chains = 1
    )
  ))

  # A posterior draws_array with one chain.
  expect_s3_class(draws, "draws_array")
  expect_equal(posterior::nchains(draws), 1)

  # The fitted parameters are the two infection probabilities plus the named
  # recovery rate gamma.
  vars <- posterior::variables(draws)
  expect_setequal(vars, c("eh_prob", "ih_prob", "gamma"))

  # Probabilities live on (0, 1) and the recovery rate must be positive.
  draws_mat <- posterior::as_draws_matrix(draws)
  expect_true(all(draws_mat[, "eh_prob"] > 0 & draws_mat[, "eh_prob"] < 1))
  expect_true(all(draws_mat[, "ih_prob"] > 0 & draws_mat[, "ih_prob"] < 1))
  expect_true(all(draws_mat[, "gamma"] > 0 & draws_mat[, "gamma"] < 1))

  # Every draw should be finite (no divergent NaNs leaking through).
  expect_true(all(is.finite(draws_mat)))
})
