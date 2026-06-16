# Integration test for data-raw/sir_res.R. SIR with a fitted recovery rate;
# always-on (no skip_on_cran) because it is the cheapest of the three fits and
# is the canonical "did you touch run_model()'s output API" guard.

sir_vars <- c("eh_prob", "ih_prob", "gamma")

test_that("sir_res bake fits and returns the expected draws_array variables", {
  draws <- fit_test_recipe(build_sir_res_spec())

  # Tier 1: the #57-class API guard. The fit must run and return a draws_array
  # whose variables are exactly the renamed SIR parameters.
  expect_s3_class(draws, "draws_array")
  expect_setequal(posterior::variables(draws), sir_vars)
})

test_that("sir_res posterior holds support invariants and matches the golden", {
  draws <- fit_test_recipe(build_sir_res_spec())

  # Tier 2 invariants: the two infection probabilities and the recovery rate
  # come back via inv_logit(), so every draw is strictly in (0, 1).
  mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(mat)))
  expect_true(all(mat > 0 & mat < 1))

  # Tier 2 regression: posterior means within a loose tolerance of the stored
  # golden (means only, backend-agnostic).
  golden <- read_golden("sir_res")
  skip_if(is.null(golden), "no stored golden for sir_res")
  means <- draws_means(draws)
  expect_setequal(names(means), names(golden))
  expect_equal(means, golden[names(means)], tolerance = golden_tol)
})
