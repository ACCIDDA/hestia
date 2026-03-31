## code to prepare `DATASET` dataset goes here
library(hestia)
library(dplyr)
source("inst/auxillary/simulation.R")

set.seed(9435)

# Simulate data from a SIIR model for 500 households with separate
# intra-household infection probabilities for each infectious compartment.
# There are three observations:
#   1. Informative about infection (e.g. culture or PCR)
#   2. second informative about recovery (e.g. IgG)
#   3. Informative about symptoms (e.g. symptom survey)
dat_sim <- sim_siir(
  eh_prob = 0.01,
  ih_prob = c(0.05, 0.025),
  n_hh = 500,
  hh_size = 1:5,
  tmax = 100,
  gamma = c(1 / 5, 1 / 3),
  split = c(0.7, 0.3),
  covs_eh = c(0, 0),
  covs_ih = c(0, 0),
  obs_prob = list(
    c(0.05, 0.95, 0.95, 0.05),
    c(0.01, 0.01, 0.01, 0.8),
    c(0.03, 1 - 1e-10, 0.03, 0.03)
  ),
  start_prob = c(1, 0, 0, 0),
  complete_enroll = TRUE
)

# Name courcome columns
siir <- dat_sim$obs |>
  rename(pcr = y1, igg = y2, symp = y3)

usethis::use_data(siir, overwrite = TRUE)
