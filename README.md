# hestia

[![R-CMD-check](https://github.com/ACCIDDA/hestia/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ACCIDDA/hestia/actions/workflows/R-CMD-check.yaml)
[![lint](https://github.com/ACCIDDA/hestia/actions/workflows/lint.yaml/badge.svg)](https://github.com/ACCIDDA/hestia/actions/workflows/lint.yaml)
[![pkgdown](https://github.com/ACCIDDA/hestia/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ACCIDDA/hestia/actions/workflows/pkgdown.yaml)
[![codecov](https://codecov.io/gh/ACCIDDA/hestia/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ACCIDDA/hestia)
[![DOI](https://zenodo.org/badge/DOI/PENDING.svg)](https://doi.org/PENDING)

<!--
Zenodo DOI badge placeholder. Zenodo tracking is already enabled for this
repo, but the DOI is only minted when the first GitHub Release is published.
After that first release, replace the two PENDING values above with the
concept DOI (the version-independent DOI that always resolves to the latest
release). Zenodo also displays the exact badge Markdown on the deposition
page once a release exists.
-->


`{hestia}` fits Bayesian compartmental infection models — such as
susceptible-infected-recovered (SIR) and multi-compartment hidden Markov
variants — from individual-level outcome data. Models are composed from
infection-process and observation-process components and fit using Stan via
`{rstan}`.

## Installation

```r
remotes::install_github("ACCIDDA/hestia")
```

Requires R >= 4.1.0 and a C++ toolchain (for Stan model compilation).

### Optional: cmdstanr backend

By default `run_model()` fits via `{rstan}` (installed automatically). It can
also fit via [`{cmdstanr}`](https://mc-stan.org/cmdstanr/), which is **not on
CRAN** and so is not installed as a dependency. To use it, install the package
and the CmdStan toolchain yourself:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
cmdstanr::install_cmdstan()
```

Then select it per call via `stan_options()`:

```r
run_model(..., stan_opts = stan_options(backend = "cmdstanr"))
```

With the cmdstanr backend you supply cmdstanr's own sampler arguments (e.g.
`parallel_chains`, `iter_warmup`); see `?stan_options`.

## Usage

Compose an infection-process model from `transmit()` (transmission between
compartments) and `progress()` (within-host progression) steps:

```r
library(hestia)

# A basic SIR model: S -> I transmission, I -> R recovery (rate fit from data)
inf_process <- make_infection_model(
  transmit(from = "S", to = "I"),
  progress(from = "I", to = "R", gamma = NA)
)
```

Pair the infection process with an observation model
(`make_observation_model()`) and fit it with `run_model()`. See the
[vignettes](https://accidda.github.io/hestia/articles/) for full worked
examples, including a multi-compartment (symptomatic/asymptomatic) model.

## Part of ACCIDDA

`{hestia}` is developed by the [Atlantic Coast Center for Infectious Disease
Dynamics and Analytics (ACCIDDA)](https://accidda.org/).
