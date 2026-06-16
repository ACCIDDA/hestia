# Integration test for data-raw/sir_cov_res.R. SIR with household covariates
# x1, x2 on both the external and internal hazards. Behind skip_on_cran() (the
# covariate model is heavier and exercises the exp() coefficient path).

cov_prob_vars <- c("eh_prob", "ih_prob", "gamma")
cov_coef_vars <- c("x1_eh", "x2_eh", "x1_ih", "x2_ih")
sir_cov_vars <- c(cov_prob_vars, cov_coef_vars)

test_that("sir_cov_res bake fits and returns the expected draws_array variables", {
  skip_on_cran()
  draws <- fit_test_recipe(build_sir_cov_res_spec(), cov = TRUE)

  # Tier 1: the #57-class API guard. The covariate path adds the four named
  # coefficients (<covcol>_eh and <covcol>_ih) to the three SIR parameters.
  expect_s3_class(draws, "draws_array")
  expect_setequal(posterior::variables(draws), sir_cov_vars)
})

test_that("sir_cov_res posterior holds support invariants and matches the golden", {
  skip_on_cran()
  draws <- fit_test_recipe(build_sir_cov_res_spec(), cov = TRUE)
  mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(mat)))

  # Tier 2 invariants, split by variable class. Probabilities and the rate come
  # back via inv_logit(): strictly in (0, 1).
  prob_mat <- mat[, cov_prob_vars, drop = FALSE]
  expect_true(all(prob_mat > 0 & prob_mat < 1))

  # Covariate coefficients come back via exp(): strictly positive and NOT
  # confined to (0, 1) (the exp() branch, not the logit one).
  for (nm in cov_coef_vars) {
    vec <- as.numeric(mat[, nm])
    expect_true(all(vec > 0), info = nm)
    expect_false(all(vec < 1), info = nm)
  }

  # Tier 2 regression: posterior means within a loose tolerance of the stored
  # golden (means only, backend-agnostic).
  golden <- read_golden("sir_cov_res")
  skip_if(is.null(golden), "no stored golden for sir_cov_res")
  means <- draws_means(draws)
  expect_setequal(names(means), names(golden))
  expect_equal(means, golden[names(means)], tolerance = golden_tol)
})
