# Tests for get_transmission_details(), which turns an infection process model
# into the matrices and parameter-indexing tables that run_model() feeds to
# Stan. These tests check the structure and the bookkeeping rather than any
# fitted values, so they run without touching Stan.

test_that("returns the documented list structure for a basic SIR model", {
  td <- get_transmission_details(sir_infection_model())

  expect_type(td, "list")
  expect_named(
    td,
    c(
      "states",
      "trans_matrix",
      "mult_matrix",
      "trans_to_fit",
      "mult_to_fit",
      "inf_states",
      "mult_inf_probs"
    )
  )
})

test_that("states are collected in order of first appearance", {
  td <- get_transmission_details(sir_infection_model())

  expect_equal(td$states, c("S", "I", "R"))
})

test_that("transition and multiplier matrices are square and labelled", {
  td <- get_transmission_details(sir_infection_model())
  n <- length(td$states)

  expect_equal(dim(td$trans_matrix), c(n, n))
  expect_equal(dim(td$mult_matrix), c(n, n))

  expect_equal(rownames(td$trans_matrix), paste("to", td$states, sep = "_"))
  expect_equal(colnames(td$trans_matrix), paste("from", td$states, sep = "_"))

  # With nothing fixed, transitions default to the tiny epsilon and the
  # split multipliers default to 1.
  expect_true(all(td$trans_matrix == 1e-10))
  expect_true(all(td$mult_matrix == 1))
})

test_that("fitted transitions are listed once per transition with a param index", {
  td <- get_transmission_details(sir_infection_model())

  # Two transitions to fit: S -> I (transmission, unnamed rate) and the named
  # gamma on I -> R.
  expect_equal(nrow(td$trans_to_fit), 2)
  expect_setequal(td$trans_to_fit$from, c("S", "I"))
  expect_setequal(td$trans_to_fit$to, c("I", "R"))

  # The only named rate (gamma) gets a positive index; the unnamed
  # transmission rate is parameterised through the infection probability and
  # so is left at 0.
  expect_equal(sum(td$trans_to_fit$param == 0), 1)
  expect_equal(sum(td$trans_to_fit$param >= 1), 1)
  gamma_row <- td$trans_to_fit[!is.na(td$trans_to_fit$rate_name), ]
  expect_equal(gamma_row$rate_name, "gamma")
  expect_true(gamma_row$param >= 1)
})

test_that("the infectious compartment is identified by index", {
  td <- get_transmission_details(sir_infection_model())

  # I is the second state, so the transmission source should index to 2.
  expect_equal(td$inf_states, which(td$states == "I"))
})

test_that("a model with no splits produces an empty multiplier table", {
  td <- get_transmission_details(sir_infection_model())

  expect_equal(nrow(td$mult_to_fit), 0)
  expect_false(isTRUE(td$mult_inf_probs))
})

test_that("a fixed transition rate is written into the matrix, not the fit table", {
  inf <- make_infection_model(
    transmit(from = "S", to = "I"),
    progress(from = "I", to = "R", gamma = 0.2)
  )
  td <- get_transmission_details(inf)

  # I -> R lands at row "to_R", column "from_I".
  expect_equal(td$trans_matrix["to_R", "from_I"], 0.2)

  # Only the unnamed transmission rate is left to fit.
  expect_equal(nrow(td$trans_to_fit), 1)
  expect_true(is.na(td$trans_to_fit$rate_name))
})

test_that("a named split is carried into the multiplier fit table", {
  td <- get_transmission_details(siir_infection_model())

  expect_equal(td$states, c("S", "Is", "Ia", "R"))

  # phi and its complementary "1-phi" partner both appear as rows to fit.
  expect_true(nrow(td$mult_to_fit) >= 1)
  expect_true("phi" %in% td$mult_to_fit$mult_name)

  # phi gets a positive parameter index.
  phi_row <- td$mult_to_fit[td$mult_to_fit$mult_name == "phi", ]
  expect_true(all(unlist(phi_row$param) >= 1))
})

test_that("separate infection probabilities expose both infectious states", {
  shared <- get_transmission_details(siir_infection_model(mult_inf_probs = FALSE))
  separate <- get_transmission_details(siir_infection_model(mult_inf_probs = TRUE))

  expect_false(isTRUE(shared$mult_inf_probs))
  expect_true(isTRUE(separate$mult_inf_probs))

  # Both Is and Ia are infectious sources regardless of the probability sharing.
  expect_setequal(separate$inf_states, which(separate$states %in% c("Is", "Ia")))
})
