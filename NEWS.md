# hestia 0.1.0

* Stan fitting and draw extraction now use `flexstanr`. Configure automatic,
  scheduler-aware within-chain threading with
  `stan_options(threading = TRUE)`.
* First release.
* Fit Bayesian compartmental infection models from individual-level outcome
  data, using Stan via 'rstan'. Models are composed from an infection process
  and an observation process and fit with a single call.
* Infection-process building blocks: `transmit()` and `progress()` define
  transitions, assembled with `make_infection_model()`. Susceptible-infected-
  recovered (SIR), split symptomatic/asymptomatic (SIIR), and similar
  multi-compartment structures are supported.
* `make_observation_model()` defines the observation process (per-compartment
  detection probabilities for one or more tests).
* `run_model()` fits the assembled model to the data and returns posterior
  draws on the natural (model) scale.
* `get_transmission_details()` exposes the parameter structure of an infection
  model.
* Bundled simulated datasets and fitted-draw fixtures for examples and the
  vignettes: `sir`, `sir_cov`, `siir`, and the `sir_res`, `sir_cov_res`,
  `siir_res` posterior draws.
* Vignettes: "SIR" and "multiple infection compartments".
