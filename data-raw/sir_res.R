## code to prepare `DATASET` dataset goes here
library(hestia)
library(dplyr)

# Basic SIR, fit recovery rate
inf_process <- make_infection_model(
  transmit(from = "S", to = "I"),
  progress(from = "I", to = "R", gamma = NA)
)

# Observation process
obs_process <- make_observation_model(
  pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
  igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
)

# Run stan model on SIR package data
sir_res <- run_model(
  inf_model = inf_process,
  obs_model = obs_process,
  data = sir,
  init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10),
  iter = 1000,
  cores = 4
)

usethis::use_data(sir_res, overwrite = TRUE)
