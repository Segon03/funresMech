# funresMech: Mechanistic Functional Response Analysis with Simulation-Based Likelihood

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/funresMech)](https://CRAN.R-project.org/package=funresMech)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> **Development version (1.1.0.9000).** This branch contains the development
> version of the forthcoming release 1.2.0. It has not been submitted to CRAN
> or archived in Zenodo; the current release on CRAN is 1.1.0.

**funresMech** fits mechanistic functional response models to experimental
counts of hosts parasitized or prey consumed. Instead of assuming a binomial or
beta-binomial distribution for the counts, it obtains their probability
distribution by simulating the foraging process (search, encounter and
handling) and estimates the parameters by maximum likelihood (Okuyama 2012,
2026). The density-scaling exponent *z*, which distinguishes type II (*z* = 1)
from type III (*z* > 1) responses, is estimated by likelihood profiling.

The package is intended for ecologists who are not specialists in
simulation-based inference. Every analysis can be run from an interactive
application, and every result comes with diagnostics that say whether it can
be trusted and what to do if it cannot.

## What the package does

* **Model specification** (`funres_model()`): distributions of search times
  (Gamma, exponential, lognormal, Weibull) and of handling times (lognormal,
  Gamma, Weibull, constant), or distributions written by the user in R or C++
  (`funres_sampler()`); parasitoid process (attacked hosts remain available)
  or predator process (attacked prey are removed); parameters held fixed or
  estimated.
* **Estimation** (`fit_funres()`, `funres_control()`): differential evolution
  with common random numbers and a Nelder-Mead polish; likelihood profile of
  *z* with a 95% likelihood-ratio interval; automatic widening of search
  bounds; comparison with *z* = 1 and with a beta-binomial model of the same
  mean by AIC. Time is handled internally in units of the trial duration, so
  results do not depend on whether *T* is given in minutes, hours or days.
* **Reliability** (`reliability()`, `verify_fit()`): checks of Monte Carlo
  noise, convergence, search bounds, the shape of the profile and the
  agreement between independent runs, each with a recommended action.
* **Data screening** (`screen_outliers()`): flags atypical trials before the
  fit, without removing them.
* **Simulation** (`simulate_response()`): simulated experiments under any
  model, for checking fits and designing experiments.
* **Interactive application** (`run_app()`): data upload, model
  specification, fitting, reliability checks, plots and an HTML report,
  without writing code.

## Installation

```r
install.packages("funresMech")                                   # release on CRAN (1.1.0)
pak::pkg_install("Segon03/funresMech@dev-frMech-1.1.0.9000")     # development version
```

Installing from source requires a C++ compiler (Rtools on Windows, Xcode
command line tools on macOS).

## Quick start

```r
library(funresMech)

# One row per trial: host density and number of hosts parasitized
d <- read.csv("my_experiment.csv")       # columns dens and par

screen_outliers(d, T = 24)               # check for atypical trials first

set.seed(1)
fit <- fit_funres(d, T = 24, time_unit = "h",
                  model = funres_model(),             # Gamma search, lognormal handling
                  control = funres_control())     # documented defaults
fit                                      # estimates, interval of z, AIC, reliability
fit <- verify_fit(fit)                   # independent runs, more effort, wider bounds
reliability(fit)                         # report only when no check fails
```

Or, without writing code:

```r
funresMech::run_app()
```

## Choosing the computational effort

The likelihood is approximated by simulation, so fits take minutes to hours.
All settings of `funres_control()` have documented defaults and can be
changed: simulated trials per density (`n_sim`, default 5,000), generations
(`itermax`, 100) and population (`NP`, 40) of differential evolution,
re-evaluations at each optimum (`n_reeval`, 5), search bounds and others.
The effort a dataset needs depends on the data and the model; the
reliability checks say whether it was enough.

As a guide, with about 60 trials at 8 to 13 densities and the profile
evaluated in parallel on 8 processor cores, a fit with `n_sim = 2000` and
`itermax = 50` takes tens of minutes, a fit with the defaults about one to
two hours, and a fit with `n_sim = 10000` and `itermax = 200` several hours.
Use `future::plan(future::multisession)` before `fit_funres()` to evaluate the
profile in parallel.

## How reliable are the results?

The estimation procedure was validated by simulation. In 288 fits to 216
datasets with known parameters (three values of *z*, three levels of search
variability and designs with 5, 8 and 13 host densities), the estimates of *z*
and *h* were unbiased, and the 95% interval of *z* contained the true value in
91% of the fits (95% confidence interval 87% to 94%). The interval is
therefore slightly narrower than nominal, mostly when search times are highly
variable, and its coverage reached 94% with 13 host densities. With the
default settings, the checks in `reliability()` gave false alarms in at most 6%
of the fits. Further simulations examined misspecified models, optimizers and
computational budgets. See the
vignette *How reliable is my fit?* and the `validation/` folder of the
repository for the scripts and results.

## Documentation

* `vignette("funresMech")`: getting started.
* `vignette("reliability", package = "funresMech")`: how reliable is my fit?
* `vignette("process-and-distributions", package = "funresMech")`: parasitoid
  or predator process, and the distributions of search and handling times.
* `vignette("units-and-bounds", package = "funresMech")`: time units, search
  bounds and open intervals.

## Citation

```r
citation("funresMech")
```

Please cite the package and the method:

* Núñez-Campero S. funresMech: mechanistic functional response analysis
  with simulation-based likelihood. R package version 1.1.0.9000
  (development version of 1.2.0).
  doi: 10.32614/CRAN.package.funresMech
* Okuyama T (2012). A likelihood approach for functional response models.
  *Biological Control* 60: 103–107. doi: 10.1016/j.biocontrol.2011.10.008
* Okuyama T (2026). Parametric assumptions in parasitoid functional response
  analysis. *Journal of Applied Entomology*. doi: 10.1111/jen.70148

## Support

Questions, bug reports and feature requests:
<https://github.com/Segon03/funresMech/issues>.

## License

MIT © Segundo Núñez-Campero
