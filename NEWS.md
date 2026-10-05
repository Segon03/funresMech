# funresMech 1.1.0

This version replaces the simulation and optimisation core. The model is the
same (Okuyama 2012, 2026); what changes is how it is computed and how the
results are summarised. It was validated against Okuyama (2026, Table S1)
with the datasets D2 and D3 of FoRAGE, and with simulated data
(parameter recovery).

## New simulation engine

* The stochastic foraging process now runs in C++ through 'Rcpp'
  (`src/motor.cpp`); with the same seed it reproduces exactly the pure-R
  reference implementation and the simulator of Okuyama (2026).
  It is about 40-65 times faster.
* Each trial is stopped once all hosts have been parasitised. The response
  is the number of distinct parasitised hosts, so this does not change the
  simulated distribution, but it removes pathological run times when the
  handling time and the Gamma shape approach zero.
* No random seed is set inside the package: `set.seed()` of the user
  controls every result, including the parallel profile.

## Changes that alter results

* **`s` is now the standard deviation of the handling time on the natural
  scale** (same units as `h`), as in Okuyama (2026), instead of the standard
  deviation on the log scale. The Lognormal still has mean `h`; the
  equivalent log-scale value is `sqrt(log(1 + (s / h)^2))`. Values of `s`
  from version 1.0.4 are not comparable.
* **The 95% interval of z uses the likelihood-ratio threshold
  `NLL_min + qchisq(0.95, 1) / 2` (+1.92).** Version 1.0.4 used +3.84, which
  gives an interval of about 99.5%.
* **Optimisation on the log scale.** The search is carried out on
  `log(a * H_ref^z)`, `log(h)`, `log(k)` and `sqrt(s)`, with wide bounds
  (1.0.4 used fixed bounds that excluded the real values of the validation
  datasets: h <= 0.5, k >= 0.5, z <= 3). A warning is shown when a parameter
  ends on a bound.
* **`z_hat` is the minimum of the likelihood profile**, and the AIC of the
  model with free z uses its NLL (previously a separate five-parameter fit,
  which could be worse than the profile because of Monte Carlo noise).
* The default z grid is 0.5-5 in steps of 0.25 (it was 0.5-1.5).

## New features

* The profile is computed in parallel over the grid of z (`future`); the
  previous `plan(multisession)` had no effect on the fit.
* If the upper limit of the interval of z is open, the grid is extended
  (z = 8, 10, 15, 20) and the limit is reported as "not identified" when the
  profile stays flat (typical when `k` is very small). When a limit falls in a
  bracket wider than the grid step, extra values of z are evaluated inside it
  so that the interpolated limit is as precise as the rest of the grid.
* Diagnostics tab and report section: small `k`, parameters on a bound,
  open intervals and noisy likelihood.
* Data screening for atypical trials before the fit (leave-one-out
  Beta-Binomial, Bonferroni threshold). Data are never removed: flagged
  trials should be checked and the analysis repeated without them.
* If the parallel workers cannot load this version of the package (for
  example, an older version is installed), the profile is evaluated
  sequentially with a message instead of failing.
* The NLL is re-evaluated at each optimum with fresh random numbers and
  reported alongside the optimiser's value (`nll_reeval`).

## Internal

* New dependency: 'Rcpp' (`Imports` and `LinkingTo`).
* `fit_full()` was replaced by `fit_profile()` and `fit_z_fixed()`.
* Package help page (`?funresMech`) and extended `?run_app`.
* New tests: engine against the R reference, saturation cut-off,
  reparameterisation, likelihood-ratio interval, profile fit and screening.
