test_that("can't multiple have observations with same name", {
  expect_error(
    obs_process <- make_observation_model(
      pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
      pcr = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
    )
  )
})

test_that("length of output list equals number of observations", {
  expect_length(
    make_observation_model(
      a = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
      b = c("S" = 0.01, "I" = 0.01, "R" = 0.8),
      c = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
    ),
    3
  )
})

test_that("entry validation rejects bad observation specifications", {
  # No observations supplied
  expect_error(
    make_observation_model(),
    "observation"
  )

  # Unnamed observation specification
  expect_error(
    make_observation_model(c("S" = 0.05, "I" = 0.95, "R" = 0.05)),
    "named"
  )

  # Probability outside [0, 1]
  expect_error(
    make_observation_model(pcr = c("S" = 0.05, "I" = 1.5, "R" = 0.05)),
    "pcr"
  )

  # Detection probabilities must themselves be named by compartment
  expect_error(
    make_observation_model(pcr = c(0.05, 0.95, 0.05)),
    "pcr"
  )
})

test_that("estimate_test validates its arguments", {
  pcr <- estimate_test(positive = "I", se = c(9, 1), sp = c(19, 1))
  expect_s3_class(pcr, "hestia_estimated_test")
  expect_identical(make_observation_model(pcr = pcr)$pcr, pcr)

  expect_error(estimate_test(positive = character(0), se = c(9, 1), sp = c(19, 1)))
  expect_error(estimate_test(positive = "I", se = c(9, 0), sp = c(19, 1)), "se")
  expect_error(estimate_test(positive = "I", se = c(9, 1), sp = 19), "sp")
})

test_that("stan data matches compartments by name", {
  inf_model <- sir_infection_model()
  data <- sir_subset(2)
  init_probs <- c(1 - 2 * 1e-10, 1e-10, 1e-10)

  ordered <- make_stan_data(
    inf_model,
    sir_observation_model(),
    data,
    init_probs
  )
  reordered <- make_stan_data(
    inf_model,
    make_observation_model(
      pcr = c("R" = 0.05, "S" = 0.05, "I" = 0.95),
      igg = c("I" = 0.01, "R" = 0.8, "S" = 0.01)
    ),
    data,
    init_probs
  )
  expect_equal(reordered$obs_fixed, ordered$obs_fixed)
  expect_equal(ordered$obs_fixed[1, 2, ], c(0.05, 0.95, 0.05))

  # Compartment names must match the infection model
  expect_error(
    make_stan_data(
      inf_model,
      make_observation_model(pcr = c("S" = 0.05, "I" = 0.95, "X" = 0.05)),
      data,
      init_probs
    ),
    "each compartment"
  )
})

test_that("stan data sets up estimated sensitivity and specificity", {
  inf_model <- sir_infection_model()
  data <- sir_subset(2)
  init_probs <- c(1 - 2 * 1e-10, 1e-10, 1e-10)

  dat_stan <- make_stan_data(
    inf_model,
    make_observation_model(
      pcr = estimate_test(positive = "I", se = c(9, 1), sp = c(19, 1)),
      igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
    ),
    data,
    init_probs
  )

  expect_equal(dat_stan$n_est, 1)
  expect_equal(as.vector(dat_stan$est_test), 1)
  expect_equal(dat_stan$pos_state, matrix(c(0, 1, 0), nrow = 1))
  expect_equal(unname(dat_stan$se_prior[1, ]), c(9, 1))
  expect_equal(unname(dat_stan$sp_prior[1, ]), c(19, 1))
  expect_equal(dat_stan$obs_fixed[2, 2, ], c(0.01, 0.01, 0.8))

  # Inits start at the prior means, inside the sens > 1 - spec region
  init <- init_obs_params(dat_stan)
  expect_equal(as.vector(init$spec), 0.95)
  expect_equal(as.vector(init$sens), 0.9)

  # Positive compartments must be some, but not all, of the compartments
  for (positive in list("X", c("S", "I", "R"))) {
    expect_error(
      make_stan_data(
        inf_model,
        make_observation_model(
          pcr = estimate_test(positive = positive, se = c(9, 1), sp = c(19, 1))
        ),
        data,
        init_probs
      ),
      "not all"
    )
  }
})
