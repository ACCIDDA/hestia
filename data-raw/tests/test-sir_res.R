# Integration test for data-raw/sir_res.R. SIR with a fitted recovery rate;
# always-on (the cheapest of the three fits) and the canonical "did you touch
# run_model()'s output API" guard. One fit covers both tiers.

test_that("sir_res bake fits, returns the expected variables, and matches golden", {
  draws <- fit_test_recipe(build_sir_res_spec())

  # Tier 1: the #57-class API guard. The fit must run and return a draws_array
  # whose variables match the bundled sir_res fixture (derived, so a bake rename
  # propagates here automatically).
  expect_s3_class(draws, "draws_array")
  expect_setequal(posterior::variables(draws), posterior::variables(sir_res))

  # Tier 2 invariants: the infection probabilities and the recovery rate come
  # back via inv_logit(), so every draw is strictly in (0, 1).
  mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(mat)))
  expect_true(all(mat > 0 & mat < 1))

  # Tier 2 regression: posterior means within tolerance of the stored golden.
  expect_means_close(draws_means(draws), read_golden("sir_res"))
})
