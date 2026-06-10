## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

* checking installed package size ... NOTE
    installed size is ~92Mb
    sub-directories of 1Mb or more: libs
  The package compiles its Stan models at install time via 'rstan'. The large
  `libs` directory is the compiled model objects and is expected for packages
  that link to Stan / StanHeaders.

## Test environments

* GitHub Actions (R-CMD-check workflow), R CMD check with `--as-cran`:
  * ubuntu-latest: R-devel, R-release, R-oldrel
  * macos-latest: R-devel, R-release, R-oldrel
  * windows-latest: R-devel, R-release, R-oldrel

<!-- TODO before submission: also run win-builder (devtools::check_win_devel()
     / check_win_release()) and R-hub (rhub::rhub_check()), and record results
     here. -->

## Downstream dependencies

There are no downstream dependencies (new package).
