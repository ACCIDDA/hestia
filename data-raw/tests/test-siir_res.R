# Integration test for data-raw/siir_res.R. SIIR with mult_inf_probs = TRUE, so
# the two infectious compartments (Is, Ia) each get their own intra-household
# infection probability. Behind skip_on_cran() (heavier than the SIR fit).

siir_vars <- c("eh_prob", "ih_prob_Is", "ih_prob_Ia", "gamma_s", "gamma_a", "phi")

test_that("siir_res bake fits and returns the expected draws_array variables", {
  skip_on_cran()
  draws <- fit_test_recipe(build_siir_res_spec())

  # Tier 1: the #57-class API guard. mult_inf_probs = TRUE must yield TWO
  # ih_prob_* variables (one per infectious state), plus both recovery rates and
  # the split phi.
  expect_s3_class(draws, "draws_array")
  expect_setequal(posterior::variables(draws), siir_vars)
})

test_that("siir_res posterior holds support invariants and matches the golden", {
  skip_on_cran()
  draws <- fit_test_recipe(build_siir_res_spec())

  # Tier 2 invariants: every SIIR parameter (probabilities, rates, split) comes
  # back via inv_logit(), so every draw is strictly in (0, 1).
  mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(mat)))
  expect_true(all(mat > 0 & mat < 1))

  # Tier 2 regression: posterior means within a loose tolerance of the stored
  # golden (means only, backend-agnostic).
  golden <- read_golden("siir_res")
  skip_if(is.null(golden), "no stored golden for siir_res")
  means <- draws_means(draws)
  expect_setequal(names(means), names(golden))
  expect_equal(means, golden[names(means)], tolerance = golden_tol)
})
