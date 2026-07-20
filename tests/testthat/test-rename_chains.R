# Minimal live covariate fit to cover the covariate post-processing path in
# rename_chains() (R/hestia_functions.R ~1062-1093): the coefficient
# subsetting/renaming and the exp() back-transform added in #19.
#
# A full covariate run_model() fit hangs the covr coverage job (-O0) for over an
# hour (see #46). This test instead fits a single short chain (chains = 1,
# iter = 100) on a tiny 3-household slice of sir_cov, mirroring imuGAP's
# test-sampling-smoke.R, so it stays fast even under covr.
#
# run_model()'s init builder hard-codes a 4-chain init list, so a chains = 1 fit
# cannot go through run_model() (it errors "initial value list mismatchs number
# of chains"). We therefore replicate run_model()'s covariate path here: build
# the stan data, supply a matching 1-chain init, sample hmm_cov directly, and
# hand the result to rename_chains() exactly as run_model() does.

test_that("rename_chains covers the covariate path on a minimal live fit", {
  skip_on_cran()

  inf_model <- sir_infection_model()
  obs_model <- sir_observation_model()

  # Tiny slice: households 1-3 and their matching covariate rows. sir_cov's
  # covariates are positional, one row per individual in (hh_id, part_id)
  # order, and households 1-3 hold the first 7 individuals.
  data_sub <- sir_cov$observations[sir_cov$observations$hh_id <= 3, ]
  cov_sub <- sir_cov$covariates[1:7, ]
  init_probs <- c(1 - 2 * 1e-10, 1e-10, 1e-10)

  dat_stan <- make_stan_data(
    inf_model, obs_model, data_sub, init_probs,
    epsilon = 1e-10, ih_cov = cov_sub, eh_cov = cov_sub
  )
  # run_model() sets these; mirror it here since we sample hmm_cov directly.
  dat_stan$save_llik <- 0L
  dat_stan$save_states <- 0L

  # One-chain init mirroring run_model()'s covariate inner list, sized from the
  # stan data so the dimensions match. Reuses the package's internal logit().
  init <- list(list(
    logit_params = array(rep(logit(0.5), dat_stan$n_params)),
    logit_mult_params = array(rep(logit(0.5), dat_stan$n_mult_params)),
    beta_eh = rep(0, dat_stan$k_eh),
    beta_ih = rep(0, dat_stan$k_ih),
    beta0_eh = logit(0.02),
    beta0_ih = array(rep(logit(0.02), dat_stan$n_inf_prob))
  ))

  fit <- suppressWarnings(rstan::sampling(
    stanmodels$hmm_cov,
    data = dat_stan,
    iter = 100,
    chains = 1,
    init = init,
    refresh = 0,
    seed = 1L
  ))

  stan_out <- list(
    stan_fit = fit,
    stan_data = dat_stan,
    eh_cov_names = colnames(cov_sub),
    ih_cov_names = colnames(cov_sub)
  )

  out <- rename_chains(inf_model, stan_out)

  expect_s3_class(out, "draws_array")
  expect_setequal(
    posterior::variables(out),
    c("eh_prob", "ih_prob", "gamma", "x1_eh", "x2_eh", "x1_ih", "x2_ih")
  )

  # Probabilities and the recovery rate come back via inv_logit(): bounded to
  # (0, 1).
  for (nm in c("eh_prob", "ih_prob", "gamma")) {
    vec <- as.numeric(out[, , nm])
    expect_true(all(vec > 0 & vec < 1), info = nm)
  }

  # Covariate coefficients come back via exp(): strictly positive and not
  # confined to (0, 1) (the exp() branch, not the logit one).
  for (nm in c("x1_eh", "x2_eh", "x1_ih", "x2_ih")) {
    vec <- as.numeric(out[, , nm])
    expect_true(all(vec > 0), info = nm)
    expect_false(all(vec < 1), info = nm)
  }
})
