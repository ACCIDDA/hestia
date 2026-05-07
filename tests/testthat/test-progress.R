test_that("one row per destination compartment", {
  expect_shape(
    progress(
      from = c("A"),
      to = c("B1", "B2", "B3"),
      split = c("boo", "baa"),
      rate1 = NA
    ),
    nrow = 3
  )
})

test_that("error when inconsistent number of to compartments and splits", {
  # Provide split when only one destination compartment
  expect_error(progress(from = "A", to = "B", split = "blah", rate = NA))

  # Don't provide split when only one destination compartment
  expect_error(progress(from = "A", to = c("B1", "B2"), rate = NA))

  # Provide too few split names
  expect_error(progress(
    from = c("A"),
    to = c("B1", "B2", "B3"),
    split = "blah",
    rate1 = NA
  ))

  # Provide too many split names
  expect_error(progress(
    from = "A",
    to = c("B1", "B2"),
    split = c("split1", "split2"),
    rate = NA
  ))

  # Provide too few split values
  expect_error(progress(
    from = c("A"),
    to = c("B1", "B2", "B3"),
    split = 0.7,
    rate1 = NA
  ))

  # Provide too many split values
  expect_error(progress(
    from = "A",
    to = c("B1", "B2"),
    split = c(0.7, 0.1),
    rate = NA
  ))

  # Split values sum to greater than 1
  expect_error(progress(
    from = c("A"),
    to = c("B1", "B2", "B3"),
    split = c(0.7, 0.5),
    rate1 = NA
  ))
})

test_that("error when no rate information provided", {
  expect_error(progress(from = "A", to = "B"))
})
