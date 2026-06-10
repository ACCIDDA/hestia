# Changelog

## hestia 0.1.0

- First release.
- Fit Bayesian compartmental infection models from individual-level
  outcome data, using Stan via ‘rstan’. Models are composed from an
  infection process and an observation process and fit with a single
  call.
- Infection-process building blocks:
  [`transmit()`](https://accidda.github.io/hestia/reference/transmit.md)
  and
  [`progress()`](https://accidda.github.io/hestia/reference/progress.md)
  define transitions, assembled with
  [`make_infection_model()`](https://accidda.github.io/hestia/reference/make_infection_model.md).
  Susceptible-infected- recovered (SIR), split symptomatic/asymptomatic
  (SIIR), and similar multi-compartment structures are supported.
- [`make_observation_model()`](https://accidda.github.io/hestia/reference/make_observation_model.md)
  defines the observation process (per-compartment detection
  probabilities for one or more tests).
- [`run_model()`](https://accidda.github.io/hestia/reference/run_model.md)
  fits the assembled model to the data and returns posterior draws on
  the natural (model) scale.
- [`get_transmission_details()`](https://accidda.github.io/hestia/reference/get_transmission_details.md)
  exposes the parameter structure of an infection model.
- Bundled simulated datasets and fitted-draw fixtures for examples and
  the vignettes: `sir`, `sir_cov`, `siir`, and the `sir_res`,
  `sir_cov_res`, `siir_res` posterior draws.
- Vignettes: “SIR” and “multiple infection compartments”.
