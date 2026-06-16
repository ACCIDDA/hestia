# Shared bake specifications for the data-raw/*_res.R scripts.
#
# Each build_*_res_spec() returns the run_model() arguments that define a fitted
# reference dataset EXCEPT the sampler budget (stan_opts). The production bake
# scripts source this file and call run_model() with the production budget; the
# integration tests under data-raw/tests/ source the same file and call
# run_model() with a tiny seeded budget. Pinning the model definition in one
# place keeps the two budgets from drifting, so a change to the API (the
# #57-class break) fails both the bake and the tests.
#
# The model objects are built with the package's own constructors, so this file
# must be sourced after the package is attached or loaded (library(hestia) in the
# bake scripts, pkgload::load_all() in the tests).

# SIR with a fitted recovery rate, observed by a PCR-like and an IgG-like test.
# Backs data-raw/sir_res.R (bundled `sir`).
build_sir_res_spec <- function() {
  list(
    inf_model = make_infection_model(
      transmit(from = "S", to = "I"),
      progress(from = "I", to = "R", gamma = NA)
    ),
    obs_model = make_observation_model(
      pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
      igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
    ),
    data = sir,
    init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10)
  )
}

# SIIR splitting infection into symptomatic (Is) and asymptomatic (Ia) via a
# fitted split phi, with separate intra-household infection probabilities per
# infectious compartment (mult_inf_probs = TRUE). Backs data-raw/siir_res.R
# (bundled `siir`).
build_siir_res_spec <- function() {
  list(
    inf_model = make_infection_model(
      transmit(
        from = "S",
        to = c("Is", "Ia"),
        source = c("Is", "Ia"),
        split = "phi"
      ),
      progress(from = "Is", to = "R", gamma_s = NA),
      progress(from = "Ia", to = "R", gamma_a = NA),
      mult_inf_probs = TRUE
    ),
    obs_model = make_observation_model(
      pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.95, "R" = 0.05),
      igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8),
      symp = c("S" = 0.03, "Is" = 1 - 1e-10, "Ia" = 0.03, "R" = 0.03)
    ),
    data = siir,
    init_probs = c(1 - 3 * 1e-10, 1e-10, 1e-10, 1e-10)
  )
}

# SIR as above but with household-level covariates x1, x2 on both the external
# and internal hazards. The same covariate frame feeds ih_cov and eh_cov, as
# make_stan_data() requires both. Backs data-raw/sir_cov_res.R (bundled
# `sir_cov`, a list of $observations and $covariates).
build_sir_cov_res_spec <- function() {
  list(
    inf_model = make_infection_model(
      transmit(from = "S", to = "I"),
      progress(from = "I", to = "R", gamma = NA)
    ),
    obs_model = make_observation_model(
      pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
      igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
    ),
    data = sir_cov$observations,
    ih_cov = sir_cov$covariates,
    eh_cov = sir_cov$covariates,
    init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10)
  )
}
