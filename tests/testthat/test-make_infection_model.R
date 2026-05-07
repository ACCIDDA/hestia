test_that("Source compartments are valid", {
  # Get an error when you define a source compartment that isn't defined
  # elsewhere in the model specification
  expect_error(make_infection_model(
    transmit(from = "S", to = "E", source = "A"),
    progress(from = "E", to = "I", sigma = NA),
    progress(from = "I", to = "R", gamma = NA)
  ))
})

test_that("No discontinuities in model structure", {
  # Get error when missing E>I transition
  expect_error(make_infection_model(
    transmit(from = "S", to = "E", source = "I"),
    progress(from = "I", to = "R", gamma = NA)
  ))
})

test_that("Requires at least one tranmit component", {
  expect_error(make_infection_model(
    progress(from = "E", to = "I", sigma = NA),
    progress(from = "I", to = "R", gamma = NA)
  ))
})
