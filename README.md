# hestia

[![R-CMD-check](https://github.com/ACCIDDA/hestia/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ACCIDDA/hestia/actions/workflows/R-CMD-check.yaml)
[![lint](https://github.com/ACCIDDA/hestia/actions/workflows/lint.yaml/badge.svg)](https://github.com/ACCIDDA/hestia/actions/workflows/lint.yaml)
[![pkgdown](https://github.com/ACCIDDA/hestia/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ACCIDDA/hestia/actions/workflows/pkgdown.yaml)
[![codecov](https://codecov.io/gh/ACCIDDA/hestia/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ACCIDDA/hestia)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20837539.svg)](https://doi.org/10.5281/zenodo.20837539)


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