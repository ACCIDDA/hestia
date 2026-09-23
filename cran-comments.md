## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission.

* checking installed package size ... NOTE
    installed size is ~92Mb
    sub-directories of 1Mb or more: libs
  The package compiles its Stan models at install time via 'rstan'. The large
  `libs` directory is the compiled model objects and is expected for packages
  that link to Stan / StanHeaders.

* checking compilation flags used ... NOTE
    Compilation used the following non-portable flag(s):
      '-Wa,-mbig-obj'
  This flag is set in `src/Makevars.win` only, so it affects Windows builds
  exclusively. The Stan-generated model code in this package produces object
  files with more COFF sections than the mingw-w64 toolchain's default limit
  (~65k) allows; without `-mbig-obj` the Windows build fails to link with a
  "file too big" error. This is a known, unavoidable limitation of the
  mingw-w64 toolchain used by R on Windows for packages that compile Stan
  models, and the same NOTE is accepted for other Stan-based packages on
  CRAN (e.g. `rstanarm`, `brms`).

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
