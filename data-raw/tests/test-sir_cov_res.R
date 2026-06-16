# Integration test for data-raw/sir_cov_res.R. SIR with household covariates on
# both the external and internal hazards. Gated behind HESTIA_RUN_INTEGRATION
# (the covariate model is heavier and exercises the exp() coefficient path).

test_that("sir_cov_res bake fits, returns the expected variables, and matches golden", {
  skip_if_not(
    run_full_integration(),
    "set HESTIA_RUN_INTEGRATION=true to run the heavier integration fits"
  )
  draws <- fit_test_recipe(build_sir_cov_res_spec(), cov = TRUE)

  # Variable classes derived from the covariate columns: one coefficient per
  # column on each hazard (exp scale); the rest are inv_logit-scale
  # probabilities and the recovery rate.
  coef_vars <- as.vector(outer(colnames(sir_cov$covariates), c("_eh", "_ih"), paste0))
  prob_vars <- setdiff(posterior::variables(sir_cov_res), coef_vars)

  # Tier 1: the #57-class API guard, against the bundled sir_cov_res fixture.
  expect_s3_class(draws, "draws_array")
  expect_setequal(posterior::variables(draws), posterior::variables(sir_cov_res))

  mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(mat)))

  # Tier 2 invariants, split by class. Probabilities and the rate are inv_logit:
  # strictly in (0, 1). Covariate coefficients are exp(): strictly positive (the
  # exp-scale magnitude itself is pinned by the golden means below).
  expect_true(all(mat[, prob_vars] > 0 & mat[, prob_vars] < 1))
  expect_true(all(mat[, coef_vars] > 0))

  # Tier 2 regression: posterior means within tolerance of the stored golden.
  expect_means_close(draws_means(draws), read_golden("sir_cov_res"))
})
