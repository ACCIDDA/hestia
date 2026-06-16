## code to prepare `siir_res` dataset goes here
library(hestia)

# Shared bake spec: model definition lives in recipes.R so the production bake
# and the data-raw/tests/ integration tests cannot drift.
source("data-raw/recipes.R")

spec <- build_siir_res_spec()

# Run the Stan model on the bundled `siir` data with the production budget.
siir_res <- run_model(
  inf_model = spec$inf_model,
  obs_model = spec$obs_model,
  data = spec$data,
  init_probs = spec$init_probs,
  stan_opts = stan_options(iter = 1000, cores = 4)
)

usethis::use_data(siir_res, overwrite = TRUE)
