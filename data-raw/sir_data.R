## code to prepare `DATASET` dataset goes here
library(hestia)
library(dplyr)
source("inst/auxillary/simulation.R")

set.seed(6842)

# Simulate data from a basic SIR model for 500 households where everyone starts
# susceptible. There are two observations:
#   1. Informative about infection (e.g. culture or PCR)
#   2. Informative about recovery (e.g. IgG)
dat_sim <- sim_sir(
  eh_prob = 0.01,
  ih_prob = 0.05,
  n_hh = 500,
  hh_size = 1:5,
  tmax = 100,
  gamma = 1 / 5,
  obs_prob = list(c(0.05, 0.95, 0.05), c(0.01, 0.01, 0.8)),
  start_prob = c(1, 0, 0),
  complete_enroll = TRUE
)

# Name outcome columns
sir <- dat_sim$obs |>
  rename(pcr = y1, igg = y2)

usethis::use_data(sir, overwrite = TRUE)
