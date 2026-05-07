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
