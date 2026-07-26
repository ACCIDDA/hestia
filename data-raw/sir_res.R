## code to prepare `sir_res` dataset goes here
library(hestia)

# Shared bake spec: model definition lives in recipes.R so the production bake
# and the data-raw/tests/ integration tests cannot drift.
source("data-raw/recipes.R")

spec <- build_sir_res_spec()

# Run the Stan model on the bundled `sir` data with the production budget.
sir_res <- run_model(
  inf_model = spec$inf_model,
  obs_model = spec$obs_model,
  data = spec$data,
  init_probs = spec$init_probs,
  stan_opts = stan_options(iter = 1000)
)

usethis::use_data(sir_res, overwrite = TRUE)
