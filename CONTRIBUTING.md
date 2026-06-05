# Contributing to hestia

Thank you for your interest in contributing to hestia! This document
explains how to report issues, propose changes, and what to expect from
the CI pipeline when you open a pull request.

## Code of Conduct

All contributors are expected to be respectful and professional in all
interactions — issues, pull requests, code reviews, and discussions.
Constructive feedback is welcome; personal attacks and dismissive
language are not.

## Reporting Bugs

Use the [bug report
template](https://github.com/ACCIDDA/hestia/issues/new?template=bug_report.yml)
when opening a new issue. It asks for a description, a minimal reprex,
your [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) output,
and your operating system.

## Suggesting Features

Use the [feature request
template](https://github.com/ACCIDDA/hestia/issues/new?template=feature_request.yml).
Describe the use case and, if possible, how the feature fits into the
existing model workflow
([`make_infection_model()`](https://accidda.github.io/hestia/reference/make_infection_model.md)
-\>
[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md)
-\> post-processing). You can also open a blank issue if neither
template fits.

## Submitting Changes

1.  **Fork & branch.** Create a feature branch off `main` (e.g.
    `fix/split-check` or `feature/missing-outcomes`).
2.  **Make your changes.** Follow the style conventions below.
3.  **Test locally.** Run `R CMD check` and the test suite before
    pushing (see the commands below).
4.  **Open a pull request** against `main` with a clear description of
    what changed and why.

### Local Development Commands

``` sh
# Full rebuild + install (required after editing any .stan file)
R CMD INSTALL --preclean .

# Regenerate documentation from roxygen comments
Rscript -e 'devtools::document()'

# Run the test suite
Rscript -e 'devtools::test()'

# Lint the package (config in .lintr)
Rscript -e 'lintr::lint_package()'

# Build and check
R CMD build . && R CMD check hestia_*.tar.gz
```

### Style

- R code is linted with `lintr` (configuration in `.lintr`).
  `R/stanmodels.R`, `R/hestia-package.R`, `inst/auxillary/`, and
  `vignettes/` are excluded from linting because they are generated or
  standalone.
- Do not hand-edit generated files: `R/stanmodels.R`,
  `src/stanExports_*.{cc,h}`, `configure`, and `configure.win` are all
  produced by
  [`rstantools::rstan_config()`](https://mc-stan.org/rstantools/reference/rstan_config.html).
  The `man/*.Rd` help pages are regenerated from roxygen comments in CI,
  so edit the roxygen, not the `.Rd`.
- Stan model variants are composed via `#include` toggling in the
  top-level `.stan` files — prefer toggling includes over duplicating
  whole model files.

## What Happens When You Open a PR

Every pull request triggers these GitHub Actions workflows:

| Workflow | What it does |
|----|----|
| **R-CMD-check** | Runs `R CMD check --as-cran` on Ubuntu, macOS, and Windows across R release, oldrel, and devel (9 jobs total). R-devel jobs are allowed to fail without blocking the PR. Stan compilation artifacts are cached per OS and R version. |
| **lint** | Runs [`lintr::lint_package()`](https://lintr.r-lib.org/reference/lint.html). |
| **test-coverage** | Runs the test suite under `covr` and uploads coverage to Codecov. |
| **pkgdown** | Builds the documentation site. On PRs this is a build-only check (no deployment); the site is deployed to GitHub Pages on pushes to `main`, published releases, and manual workflow dispatches. |

These workflows should pass before a PR is merged.

## Questions?

If something is unclear, open an issue or comment on an existing one —
we are happy to help.
