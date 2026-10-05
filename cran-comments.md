## funresMech 1.1.0 (update of 1.0.4)

This is an update of the version on CRAN (1.0.4). It replaces the simulation
and optimisation core of the package; the statistical model is unchanged
(Okuyama 2012, 2026). Details are in NEWS.md.

* The stochastic simulator is now written in C++ with 'Rcpp' (new
  `Imports`/`LinkingTo` dependency: 'Rcpp'). The package therefore needs
  compilation (`src/motor.cpp`).
* Results can differ from 1.0.4 by design: `s` is now the standard deviation
  of the handling time on the natural scale, the 95% interval of `z` uses the
  likelihood-ratio threshold of 1.92 (1.0.4 used 3.84) and the optimisation is
  carried out on the log scale with wider bounds. These changes are
  documented in NEWS.md.
* No `set.seed()` is called inside the package; all randomness is controlled
  by the caller. The package writes nothing outside `tempdir()`.
* Tests: `tests/testthat/test-modelos.R` checks the C++ engine against a
  pure-R reference implementation, the likelihood-ratio interval, the profile
  fit and the screening of atypical trials. The slow tests use
  `skip_on_cran()`.

## Test environments

* Local: Windows 11 x64, R 4.6.1 (2026-06-24 ucrt), Rtools45 (GCC 14.3.0).
* win-builder, R-devel (R Under development, 2026-09-30 r90605 ucrt).

## R CMD check results

Local, `devtools::check(args = "--as-cran")`:

0 errors | 0 warnings | 0 notes

`devtools::test()`: 79 expectations passed, 0 failures, 0 warnings.

win-builder, R-devel (R Under development (unstable) (2026-09-30 r90605 ucrt)),
`devtools::check_win_devel()`:

Status: OK
