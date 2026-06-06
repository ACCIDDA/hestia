test_that("return one set of chains per parameter", {
  # Subset to first ten households
  siir_sub <- sir[sir$hh_id <= 10, ]

  # Infection process model
  inf_mod <- make_infection_model(
    transmit(from = "S", to = c("Is", "Ia"), split = "phi"),
    progress(from = "Is", to = "R", gamma = NA),
    progress(from = "Ia", to = "R", gamma = NA)
  )

  obs_mod <- make_observation_model(
    pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.05, "R" = 0.05),
    igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8)
  )

  suppressWarnings(suppressMessages(
    model_out <- run_model(
      inf_model = inf_mod,
      obs_model = obs_mod,
      data = siir_sub,
      init_probs = c(1 - 3 * 1e-10, 1e-10, 1e-10, 1e-10),
      iter = 10
    )
  ))

  # Variable is the third dimension in a draws array_object
  # Should be 4 variables: eh_prob, ih_prob, gamma, phi
  expect_equal(length(posterior::variables(model_out)), 4)
})

test_that("multiple infections probabilities supported", {
  # Subset to first ten households
  siir_sub <- sir[sir$hh_id <= 10, ]

  # Infection process model
  inf_mod <- make_infection_model(
    transmit(from = "S", to = c("Is", "Ia"), split = "phi"),
    progress(from = "Is", to = "R", gamma = NA),
    progress(from = "Ia", to = "R", gamma = NA),
    mult_inf_probs = TRUE
  )

  obs_mod <- make_observation_model(
    pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.05, "R" = 0.05),
    igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8)
  )

  suppressWarnings(suppressMessages(
    model_out <- run_model(
      inf_model = inf_mod,
      obs_model = obs_mod,
      data = siir_sub,
      init_probs = c(1 - 3 * 1e-10, 1e-10, 1e-10, 1e-10),
      iter = 10
    )
  ))

  # Variable is the third dimension in a draws array_object
  # Should be 5 variables: eh_prob, ih_prob_Is, ih_prob_Ia, gamma, phi
  expect_equal(length(posterior::variables(model_out)), 5)

  # Should be two variables that match ih_prob
  expect_equal(length(grep("ih_prob", posterior::variables(model_out))), 2)
})

test_that("entry validation rejects bad run_model inputs before sampling", {
  inf_mod <- make_infection_model(
    transmit(from = "S", to = "I"),
    progress(from = "I", to = "R", gamma = NA)
  )
  obs_mod <- make_observation_model(
    pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
    igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
  )
  good_data <- sir[sir$hh_id <= 10, ]
  good_init <- c(1 - 2 * 1e-10, 1e-10, 1e-10)

  # init_probs that do not sum to 1
  expect_error(
    run_model(inf_mod, obs_mod, good_data, init_probs = c(0.2, 0.2, 0.2)),
    "sum to 1"
  )

  # data missing a required structural column
  expect_error(
    run_model(
      inf_mod,
      obs_mod,
      good_data[, setdiff(names(good_data), "hh_size")],
      init_probs = good_init
    ),
    "hh_size"
  )

  # data missing an outcome column named in the observation model
  expect_error(
    run_model(
      inf_mod,
      obs_mod,
      good_data[, setdiff(names(good_data), "igg")],
      init_probs = good_init
    ),
    "igg"
  )

  # obs_model that is not a named list
  expect_error(
    run_model(inf_mod, list(), good_data, init_probs = good_init),
    "obs_model"
  )
})

# Variable names match inputs
