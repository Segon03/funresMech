# funresMech 1.1.0.9000

Development version of the forthcoming release 1.2.0, available on the GitHub
branch `dev-frMech-1.1.0.9000`. It has not been submitted to CRAN or archived
in Zenodo; the current release on CRAN is 1.1.0. The changes below will be
released as version 1.2.0.

This version turns the package into a general tool for simulation-based
likelihood analysis of functional responses. The foraging process can now be
specified by the user, every fit is accompanied by reliability checks, and the
estimation procedure was validated with simulated data under a wide range of
parameters, experimental designs and model misspecifications. The default
model (Gamma search times, lognormal handling times, parasitoid process)
reproduces version 1.1.0 exactly with the same random seed.

## Model specification

* New `funres_model()`: distributions of search times (`"gamma"`,
  `"exponential"`, `"lognormal"`, `"weibull"`) and of handling times
  (`"lognormal"`, `"gamma"`, `"weibull"`, `"constant"`), chosen
  independently. For every built-in family the mean search time is
  1/(a x^z) and its dispersion is set by k (CV^2 = 1/k); handling times
  have mean h and standard deviation s.
* New `funres_sampler()` and `test_sampler()`: distributions written by the
  user, as R functions or compiled C++ functions.
* New `process = "predator"`: attacked prey are removed. For a parasitoid,
  this represents perfect and cost-free discrimination of parasitized hosts.
  The default `process = "parasitoid"` keeps attacked hosts available.
* Parameters `k`, `s` and `z` can be held fixed (`fix = list(...)`); the
  number of free parameters is propagated to the AIC.

## Estimation

* New `fit_funres()` with `funres_control()`. Three presets: `"exploratory"`
  (not for reporting), `"standard"` (default) and `"final"`.
* Default optimizer: differential evolution with common random numbers and a
  Nelder-Mead polish, chosen in a comparison of strategies with equal budgets
  on simulated data. CMA-ES (`method = "CMAES"`, package 'cmaes') and CRS2
  (`method = "CRS2"`, package 'nloptr') are available as alternatives.
* Time is handled internally in units of the trial duration `T`, so the
  results and the default search bounds do not depend on the time unit of
  the data. Bounds can be changed with `funres_control(bounds = ...)` and are
  widened automatically when an estimate ends on one of them.
* Joint estimation of `z` with the other parameters (`z_method = "joint"`) is
  available as an experimental option; the likelihood profile remains the
  default because its intervals had better coverage in simulations.
* The AIC of the fitted model is compared with that of the model with `z = 1`
  and with a beta-binomial model with the same mean.

## Interval of z

* The profile and the interval use the mean of independent re-evaluations of
  the likelihood at each optimum, not the minimum found by the optimizer.
* When an observed outcome is never produced by the simulation (probability
  floor), the re-evaluations are repeated with four times more simulations
  (`floor_retry`); points where this persists are left out of the interval
  and reported.
* Extra profile points are added around the estimate and where the profile
  crosses the 95% threshold, and the limits are interpolated with a
  shape-preserving cubic spline.
* The grid of `z` is extended upward (8, 10, 15, 20) or downward (0.25, 0.1)
  when a limit is open, and a limit that stays open after the extension is
  reported as not identified by the data.

## Reliability

* New `reliability()`: checks of Monte Carlo noise, the gap between the
  optimizer value and the re-evaluated likelihood, convergence, search
  bounds, the regularity of the profile, open intervals, extreme search
  dispersion and the preset, each with a status and a recommended action.
* New `verify_fit()`: independent runs, a run with more effort and a run with
  wider bounds at the estimate of `z`.
* New exported `screen_outliers()`: leave-one-out screening for atypical
  trials. Data are never removed.
* New `simulate_response()`: simulated experiments under any model.

## Interactive application

* Model specification panel (distributions, user samplers with a test
  button, process, fixed parameters), time unit selector, presets with an
  estimate of the run time, random seed, and editable search bounds.
* Reliability tab with the checks and a "Run verification" button; a
  permanent notice is shown while a fit is not verified.
* Model comparison table, downloads of estimates, profiles, checks and fitted
  objects, an R script that repeats the analysis, and an HTML report.
* The application no longer refers to the method by an author's name.

## Changes that alter results

* The default preset is `"standard"` (5,000 simulations per density,
  100 generations, 5 re-evaluations).
* The default grid of `z` is 0.5 to 5 in steps of 0.25, with automatic
  extension and refinement.

## Internal

* The package title is now "Mechanistic Functional Response Analysis with
  Simulation-Based Likelihood".
* 399 tests, including the exact reproduction of version 1.1.0 for the
  default model and the invariance of complete fits to the time unit.

# funresMech 1.1.0

This version replaced the simulation and optimization core. The model is the
same (Okuyama 2012, 2026); what changed is how it is computed and how the
results are summarized. It was validated with the datasets D2 and D3 of the
FoRAGE database and with simulated data (parameter recovery).

## New simulation engine

* The stochastic foraging process runs in C++ through 'Rcpp'. With the same
  seed it reproduces exactly the pure R reference implementation and the
  supplementary code of Okuyama (2026). It is about 40 to 65 times faster.
* Each trial is stopped once all hosts have been parasitized. The response
  is the number of distinct parasitized hosts, so this does not change the
  simulated distribution, but it removes pathological run times when the
  handling time and the Gamma shape approach zero.
* No random seed is set inside the package: `set.seed()` controls every
  result, including the parallel profile.

## Changes that alter results

* **`s` is the standard deviation of the handling time on the natural
  scale** (same units as `h`), instead of the standard deviation on the log
  scale. The lognormal still has mean `h`; the equivalent log-scale value is
  `sqrt(log(1 + (s / h)^2))`. Values of `s` from version 1.0.4 are not
  comparable.
* **The 95% interval of z uses the likelihood-ratio threshold
  `NLL_min + qchisq(0.95, 1) / 2` (+1.92).** Version 1.0.4 used +3.84, which
  gives an interval of about 99.5%.
* **Optimization on the log scale**, with wide bounds. Version 1.0.4 used
  fixed bounds that excluded the estimates of the validation datasets
  (h <= 0.5, k >= 0.5, z <= 3). A warning is shown when a parameter ends on a
  bound.
* **`z_hat` is the minimum of the likelihood profile**, and the AIC of the
  model with free z uses its negative log-likelihood.
* The default z grid is 0.5 to 5 in steps of 0.25 (it was 0.5 to 1.5).

## New features

* The profile is computed in parallel over the grid of z ('future').
* Open upper limits of the interval of z trigger an extension of the grid
  (z = 8, 10, 15, 20); wide brackets around a limit are refined.
* Diagnostics tab and report section: small `k`, parameters on a bound,
  open intervals and noisy likelihood.
* Screening for atypical trials before the fit (leave-one-out
  beta-binomial, Bonferroni threshold). Data are never removed.
* The negative log-likelihood is re-evaluated at each optimum with fresh
  random numbers and reported with the optimizer's value.

## Internal

* New dependency: 'Rcpp' (`Imports` and `LinkingTo`).
* Package help page (`?funresMech`) and extended `?run_app`.
* New tests: engine against the R reference, saturation stop,
  reparameterization, likelihood-ratio interval, profile fit and screening.
