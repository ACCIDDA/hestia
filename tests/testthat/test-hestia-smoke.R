# End-to-end smoke test. This fits a real (if tiny) Stan model on the bundled
# `sir` data and checks that the pipeline returns sensible, well-named draws.
# It deliberately uses a single short chain so it is cheap enough to run in CI
# on every push (no skip_on_cran()).

test_that("the package and its bundled data are available", {
  expect_true(requireNamespace("hestia", quietly = TRUE))

  data(sir, package = "hestia", envir = environment())
  expect_s3_class(sir, "data.frame")
  expected_cols <- c("t", "part_id", "hh_size", "hh_id", "pcr", "igg")
  expect_true(all(expected_cols %in% names(sir)))
})

test_that("run_model fits an SIR model end-to-end and returns named draws", {
  inf_mod <- sir_infection_model()
  obs_mod <- sir_observation_model()
  dat <- sir_subset(10)

  # run_model's default init list is built for the default 4 chains, so we keep
  # the default chain count and just shorten the run to stay cheap in CI.
  suppressWarnings(suppressMessages(
    draws <- run_model(
      inf_model = inf_mod,
      obs_model = obs_mod,
      data = dat,
      init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10),
      iter = 200
    )
  ))

  expect_s3_class(draws, "draws_array")
  expect_equal(posterior::nchains(draws), 4)

  # rename_chains labels the fitted parameters: the two infection probabilities
  # and the recovery rate gamma. They are returned on the natural (model) scale
  # (#19), so every value is a probability in (0, 1).
  param_names <- c("eh_prob", "ih_prob", "gamma")
  expect_setequal(posterior::variables(draws), param_names)
  draws_mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(draws_mat)))
  expect_true(all(draws_mat > 0 & draws_mat < 1))
})

test_that("run_model fits a covariate model and exp-transforms coefficients", {
  # Exercises the covariate code paths in run_model() and rename_chains() (#19):
  # the covariate init, the covariate variable subsetting/renaming, and the
  # exp() coefficient transform that the non-covariate smoke test never reaches.
  inf_mod <- sir_infection_model()
  obs_mod <- sir_observation_model()

  suppressWarnings(suppressMessages(
    draws <- run_model(
      inf_model = inf_mod,
      obs_model = obs_mod,
      data = sir_cov$observations,
      ih_cov = sir_cov$covariates,
      eh_cov = sir_cov$covariates,
      init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10),
      iter = 200
    )
  ))

  expect_s3_class(draws, "draws_array")
  vars <- posterior::variables(draws)

  # Base parameters come back on the natural (probability) scale.
  base_params <- c("eh_prob", "ih_prob", "gamma")
  expect_true(all(base_params %in% vars))

  # Covariate coefficients are exp()-transformed off the log scale, so they are
  # positive but not bounded to (0, 1).
  coef_params <- c("x1_eh", "x2_eh", "x1_ih", "x2_ih")
  expect_true(all(coef_params %in% vars))

  draws_mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(draws_mat)))
  expect_true(all(draws_mat[, base_params] > 0 & draws_mat[, base_params] < 1))
  expect_true(all(draws_mat[, coef_params] > 0))
})
