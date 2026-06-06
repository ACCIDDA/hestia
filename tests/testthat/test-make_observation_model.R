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
