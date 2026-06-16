## Maintainer-only: regenerate data-raw/tests/golden-means.rds
##
## The Tier-2 regression in the data-raw integration tests compares posterior
## MEANS (no draws) against a small stored golden. Regenerate it deliberately
## when the recipes or the sampler intentionally change, never as part of CI.
##
## Run from the package root after the package compiles:
##   PKG_BUILD_EXTRA_FLAGS=false Rscript -e 'pkgload::load_all(".")' \
##     -e 'source("data-raw/regenerate-golden.R")'
## or interactively: pkgload::load_all("."); source("data-raw/regenerate-golden.R")
##
## The golden is produced with the SAME tiny seeded budget the tests use
## (test_stan_opts, ~3-household slice), so seeded determinism keeps the live
## test means within golden_tol of these stored values.

source("data-raw/tests/helper-bake.R")

golden <- list(
  sir_res = draws_means(fit_test_recipe(build_sir_res_spec())),
  siir_res = draws_means(fit_test_recipe(build_siir_res_spec())),
  sir_cov_res = draws_means(fit_test_recipe(build_sir_cov_res_spec(), cov = TRUE))
)

saveRDS(golden, data_raw_path("tests", "golden-means.rds"))

message("Wrote data-raw/tests/golden-means.rds:")
for (nm in names(golden)) {
  message("  ", nm, ": ", paste(names(golden[[nm]]), collapse = ", "))
}
