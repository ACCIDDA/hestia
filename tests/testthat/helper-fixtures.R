# Shared fixtures for the hestia test suite.
#
# These helpers build the small, reusable model objects that several tests
# need so the individual test files stay focused on behaviour rather than
# setup. Everything here is cheap to construct (no Stan, no sampling).

# A basic SIR infection process: S -> I (transmission) and I -> R at a fitted
# rate gamma. This is the smallest model that still exercises both a transmit
# and a progress transition.
sir_infection_model <- function() {
  make_infection_model(
    transmit(from = "S", to = "I"),
    progress(from = "I", to = "R", gamma = NA)
  )
}

# An SIIR infection process that splits infection into symptomatic (Is) and
# asymptomatic (Ia) compartments via a fitted split parameter phi. Useful for
# checking split / multiplier handling. Set mult_inf_probs to TRUE to give the
# two infectious compartments separate intra-household infection probabilities.
siir_infection_model <- function(mult_inf_probs = FALSE) {
  make_infection_model(
    transmit(from = "S", to = c("Is", "Ia"), split = "phi"),
    progress(from = "Is", to = "R", gamma = NA),
    progress(from = "Ia", to = "R", gamma = NA),
    mult_inf_probs = mult_inf_probs
  )
}

# Observation model matched to the SIR compartments: a PCR-like test that is
# sensitive while infectious and an IgG-like test that turns on after recovery.
sir_observation_model <- function() {
  make_observation_model(
    pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
    igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
  )
}

# Observation model matched to the SIIR compartments.
siir_observation_model <- function() {
  make_observation_model(
    pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.05, "R" = 0.05),
    igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8)
  )
}

# A small slice of the bundled `sir` dataset, kept tiny so end-to-end tests
# stay fast. `n_hh` households are taken from the front of the data.
sir_subset <- function(n_hh = 10) {
  sir[sir$hh_id <= n_hh, ]
}
