from_ref <- "A"
to_ref <- c("B1", "B2", "B3")
source_ref <- c("C1", "C2")
split_name_ref <- c("boo", "baa")
split_num_ref <- c(0.7, 0.2)

test_that("One to one, no source specified, no split specified works", {
  transmit_test <- transmit(
    from = from_ref,
    to = to_ref[1]
  )

  # Should be a data frame
  expect_s3_class(transmit_test, "data.frame")

  # Should only have one row
  expect_shape(
    transmit_test,
    nrow = length(to_ref[1])
  )

  # Source should autofill from "to" argument
  expect_equal(unlist(transmit_test$source), transmit_test$to)

  # Should have expected column names
  checkmate::expect_names(
    names(transmit_test),
    permutation.of = c("from", "to", "split_name", "split_value", "source")
  )
})

test_that("One to many, no source specified, character split works", {
  transmit_test <- transmit(
    from = from_ref,
    to = to_ref,
    split = split_name_ref
  )

  # Should be a data frame
  expect_s3_class(transmit_test, "data.frame")

  # Should only have number of rows equal to length of "to"
  expect_shape(
    transmit_test,
    nrow = length(to_ref)
  )

  # Source should autofill from "to" argument
  expect_equal(transmit_test$source, rep(list(to_ref), length(to_ref)))

  # Split names fill in correctly
  expect_equal(
    transmit_test$split_name,
    c(
      split_name_ref,
      paste0("1-", paste(split_name_ref, sep = "-", collapse = "-"))
    )
  )

  # Split values are NA
  expect_equal(sum(is.na(transmit_test$split_value)), nrow(transmit_test))

  # Should have expected column names
  checkmate::expect_names(
    names(transmit_test),
    permutation.of = c("from", "to", "split_name", "split_value", "source")
  )
})


test_that("One to many, no source specified, numeric split works", {
  transmit_test <- transmit(
    from = from_ref,
    to = to_ref,
    split = split_num_ref
  )

  # Should be a data frame
  expect_s3_class(transmit_test, "data.frame")

  # Should only have number of rows equal to length of "to"
  expect_shape(
    transmit_test,
    nrow = length(to_ref)
  )

  # Source should autofill from "to" argument
  expect_equal(transmit_test$source, rep(list(to_ref), length(to_ref)))

  # Split values fill in correctly
  expect_equal(
    transmit_test$split_value,
    c(split_num_ref, 1 - sum(split_num_ref))
  )

  # Split names are NA
  expect_equal(sum(is.na(transmit_test$split_name)), nrow(transmit_test))

  # Should have expected column names
  checkmate::expect_names(
    names(transmit_test),
    permutation.of = c("from", "to", "split_name", "split_value", "source")
  )
})

test_that("One to many, source specified, numeric split works", {
  transmit_test <- transmit(
    from = from_ref,
    to = to_ref,
    split = split_num_ref,
    source = source_ref
  )

  # Should be a data frame
  expect_s3_class(transmit_test, "data.frame")

  # Should only have number of rows equal to length of "to"
  expect_shape(
    transmit_test,
    nrow = length(to_ref)
  )

  # Source should be equal to input vector for all rows
  expect_equal(transmit_test$source, rep(list(source_ref), length(to_ref)))

  # Split values fill in correctly
  expect_equal(
    transmit_test$split_value,
    c(split_num_ref, 1 - sum(split_num_ref))
  )

  # Split names are NA
  expect_equal(sum(is.na(transmit_test$split_name)), nrow(transmit_test))

  # Should have expected column names
  checkmate::expect_names(
    names(transmit_test),
    permutation.of = c("from", "to", "split_name", "split_value", "source")
  )
})

test_that("error when inconsistencies in splits", {
  # Provide split when only one destination compartment
  expect_error(transmit(from = "A", to = "B", split = "blah"))

  # Don't provide split when more than one destination compartment
  expect_error(transmit(from = "A", to = c("B1", "B2")))

  # Provide too few split names
  expect_error(progress(
    from = from_ref,
    to = to_ref,
    split = split_name_ref[1]
  ))

  # Provide too many split names
  expect_error(transmit(
    from = from_ref,
    to = to_ref,
    split = c(split_name_ref, "blah")
  ))

  # Provide too few split values
  expect_error(transmit(
    from = from_ref,
    to = to_ref,
    split = split_num_ref[1]
  ))

  # Provide too many split values
  expect_error(transmit(
    from = from_ref,
    to = to_ref,
    split = c(split_num_ref, 0.01, 0.01)
  ))

  # Split values sum to greater than 1
  expect_error(transmit(
    from = from_ref,
    to = to_ref,
    split = split_num_ref + 0.2
  ))

  # Splits not character or numeric
  expect_error(transmit(
    from = from_ref,
    to = to_ref,
    split = rep(TRUE, length(to_ref) - 1)
  ))

  # Some splits are negative
  expect_error(transmit(
    from = from_ref,
    to = to_ref,
    split = -1 * split_num_ref
  ))
})

test_that("error when empty or improper values given to 'to' or 'from'", {
  # Empty from
  expect_error(transmit(
    from = c(),
    to = to_ref
  ))

  # Empty to
  expect_error(transmit(
    from = from_ref,
    to = c()
  ))

  # Non-character from
  expect_error(transmit(
    from = 1,
    to = to_ref
  ))

  # Non-character to
  expect_error(transmit(
    from = from_ref,
    to = 1
  ))

  # More than one value for from
  expect_error(transmit(
    from = c(from_ref, "X"),
    to = to_ref
  ))
})

test_that("to and from can't overlap", {
  expect_error(transmit(
    from = from_ref,
    to = c(from_ref[1], to_ref),
    split = rep(TRUE, length(to_ref) - 1)
  ))
})

test_that("entry validation rejects bad from/to/split types", {
  # Non-character from is caught by the entry checker
  expect_error(
    transmit(from = 1, to = to_ref),
    "from"
  )

  # Missing values in to
  expect_error(
    transmit(from = from_ref, to = c("B1", NA)),
    "to"
  )

  # split that is neither NA, character, nor numeric
  expect_error(
    transmit(from = from_ref, to = to_ref, split = list(1, 2)),
    "split"
  )

  # numeric split containing NA (beyond the single-NA default)
  expect_error(
    transmit(from = from_ref, to = to_ref, split = c(0.5, NA)),
    "split"
  )
})

test_that("transmit rejects a from compartment that is also a destination", {
  # A compartment cannot transition to itself.
  expect_error(transmit(from = "S", to = "S"), "different compartments")
  expect_error(
    transmit(from = "S", to = c("S", "I"), split = "phi"),
    "different compartments"
  )
})
