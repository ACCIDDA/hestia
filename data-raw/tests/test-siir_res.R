# Integration test for data-raw/siir_res.R. SIIR with mult_inf_probs = TRUE, so
# the two infectious compartments (Is, Ia) each get their own intra-household
# infection probability. Gated behind HESTIA_RUN_INTEGRATION (heavier than SIR).

test_that("siir_res bake fits, returns the expected variables, and matches golden", {
  skip_if_not(
    run_full_integration(),
    "set HESTIA_RUN_INTEGRATION=true to run the heavier integration fits"
  )
  draws <- fit_test_recipe(build_siir_res_spec())

  # Tier 1: mult_inf_probs = TRUE must yield two ih_prob_* variables; the full
  # set is derived from the bundled siir_res fixture.
  expect_s3_class(draws, "draws_array")
  expect_setequal(posterior::variables(draws), posterior::variables(siir_res))

  # Tier 2 invariants: every SIIR parameter comes back via inv_logit(), so every
  # draw is strictly in (0, 1).
  mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(mat)))
  expect_true(all(mat > 0 & mat < 1))

  # Tier 2 regression: posterior means within tolerance of the stored golden.
  expect_means_close(draws_means(draws), read_golden("siir_res"))
})
