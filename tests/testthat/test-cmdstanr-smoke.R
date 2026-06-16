# End-to-end smoke test for the cmdstanr backend. Skips unless cmdstanr and a
# working CmdStan toolchain are available, so it is safe on CI and CRAN (which
# don't have them) and runs locally for anyone set up for cmdstanr.

test_that("run_model fits an SIR model via the cmdstanr backend", {
  skip_if_not_installed("cmdstanr")
  skip_if(
    inherits(try(cmdstanr::cmdstan_version(), silent = TRUE), "try-error"),
    "CmdStan toolchain not available"
  )

  inf_mod <- sir_infection_model()
  obs_mod <- sir_observation_model()
  dat <- sir_subset(10)

  suppressWarnings(suppressMessages(
    draws <- run_model(
      inf_model = inf_mod,
      obs_model = obs_mod,
      data = dat,
      init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10),
      stan_opts = stan_options(
        backend = "cmdstanr",
        iter_warmup = 100,
        iter_sampling = 100
      )
    )
  ))

  expect_s3_class(draws, "draws_array")
  expect_setequal(
    posterior::variables(draws),
    c("eh_prob", "ih_prob", "gamma")
  )
  draws_mat <- posterior::as_draws_matrix(draws)
  expect_true(all(is.finite(draws_mat)))
  expect_true(all(draws_mat > 0 & draws_mat < 1))
})
