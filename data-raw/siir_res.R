## code to prepare `DATASET` dataset goes here
library(hestia)
library(dplyr)

# Infection process model
inf_process <- make_infection_model(
  transmit(
    from = "S",
    to = c("Is", "Ia"),
    source = c("Is", "Ia"),
    split = "phi"
  ),
  progress(from = "Is", to = "R", gamma_s = NA),
  progress(from = "Ia", to = "R", gamma_a = NA),
  mult_inf_probs = TRUE
)


# Observation process
obs_process <- make_observation_model(
  pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.95, "R" = 0.05),
  igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8),
  symp = c("S" = 0.03, "Is" = 1 - 1e-10, "Ia" = 0.03, "R" = 0.03)
)

# Run stan model on SIR package data
siir_res <- run_model(
  inf_model = inf_process,
  obs_model = obs_process,
  data = siir,
  init_probs = c(1 - 3 * 1e-10, 1e-10, 1e-10, 1e-10),
  iter = 1000,
  cores = 4
)

usethis::use_data(siir_res, overwrite = TRUE)
