## funresMech 1.2.0 (update of 1.1.0)

This is an update of the version on CRAN (1.1.0). Details are in NEWS.md.

* New exported functions to specify the foraging process (`funres_model()`,
  `funres_sampler()`, `test_sampler()`), fit and verify models
  (`fit_funres()`, `funres_control()`, `verify_fit()`, `reliability()`),
  screen the data (`screen_outliers()`), and simulate experiments
  (`simulate_response()`).
* Four vignettes (built with 'knitr'; the code chunks that take minutes are
  not evaluated).
* Results can differ from 1.1.0 by design: the default settings, the search
  bounds (now independent of the time unit), the optimizer (differential
  evolution with common random numbers and a Nelder-Mead polish), and the
  construction of the interval of `z` changed. These changes are listed in
  NEWS.md under "Changes that alter results". The default model reproduces
  the simulation engine of 1.1.0 exactly with the same seed.
* New suggested packages: 'cmaes' and 'nloptr' (alternative optimizers, used
  only when requested) and 'knitr' (vignettes).
* No `set.seed()` is called inside the package; all randomness is controlled
  by the caller. The package writes nothing outside `tempdir()`.
* The slow tests use `skip_on_cran()`.

## Test environments

* Linux (Ubuntu), R 4.3.3: `R CMD check --as-cran`.
* To be completed before submission: local Windows 11 x64, R 4.6.1, Rtools45;
  win-builder, R-devel.

## R CMD check results

Linux, R 4.3.3: 0 errors | 0 warnings. The notes were due to the check
environment (no internet access for URL checks, no system clock
verification).

To be completed before submission with the Windows and win-builder results.
