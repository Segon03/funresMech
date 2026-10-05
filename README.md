# funresMech: Mechanistic Functional Response Analysis using the Okuyama Model

[![CRAN status](https://www.r-pkg.org/badges/version/funresMech)](https://CRAN.R-project.org/package=funresMech)
[![R-CMD-check](https://github.com/Segon03/funresMech/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Segon03/funresMech/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**funresMech** implements the mechanistic, stochastic functional response model proposed by Okuyama (2012) and extended to parasitoids in Okuyama (2026). Unlike traditional approaches that rely on heuristic distributions (binomial or beta-binomial), this package simulates the underlying search-encounter-handling process to generate the probability distribution of the data, providing a more flexible and mechanistically sound framework for functional response analysis.

The package includes:
- An interactive **Shiny application** for data exploration and model fitting.
- **Maximum likelihood estimation** using a simulation-based likelihood.
- **Likelihood profiles** for the density-scaling exponent *z*.
- **Model comparison** via AIC between full (*z* free) and restricted (\\(z = 1\\)) models.
- **Comprehensive diagnostic plots** including stochastic curves, histograms, density plots, boxplots, violins, and fan plots.

## Installation

### From CRAN
```r
install.packages("funresMech")
```

### From GitHub (development version)
```r
# Using pak (recommended)
install.packages("pak")
pak::pkg_install("Segon03/funresMech")

# Or using devtools
install.packages("devtools")
devtools::install_github("Segon03/funresMech")
```

Since version 1.1.0 the package contains compiled code (C++ via 'Rcpp'). On
Windows and macOS, installing from GitHub requires a compiler (Rtools /
Xcode command line tools); the CRAN version is distributed as a binary.

## Usage

### Launch the Shiny app
```r
library(funresMech)
run_app()
```

In the app you can:

1. Upload your dataset (CSV): one row per trial with species, host density and number of parasitised hosts.
2. Select the columns, the experiment duration *T* and the grid of *z*.
3. Run the analysis: parameters, likelihood profile of *z* with its 95% interval, AIC (free *z* vs *z* = 1), stochastic curves and distributions.
4. Review the **Diagnostics** and **Data screening** tabs, and download an HTML report.

For reproducible results call `set.seed()` before `run_app()`; the likelihood profile is computed in parallel and gives the same result as a sequential run with the same seed.

## Features

- **Mechanistic simulation:** search times follow a Gamma distribution; handling times follow a Lognormal distribution with mean `h` and standard deviation `s` (natural scale, as in Okuyama 2026). The engine is written in C++ ('Rcpp').
- **Simulated likelihood:** the probability distribution of parasitism is generated through repeated simulations.
- **Flexible density scaling:** the exponent *z* allows Type I, II and III-like responses.
- **Uncertainty:** 95% interval of *z* from the likelihood profile (likelihood-ratio threshold of 1.92 log-likelihood units), with automatic extension and refinement of the grid when a limit is open.
- **Diagnostics and screening of atypical trials** (data are never removed).
- **Interactive visualisation** with 'plotly' and **reports** in HTML.

## What is new in 1.1.0

Version 1.1.0 replaces the computational core and changes some results with
respect to 1.0.4 (scale of `s`, threshold of the interval, search on the log
scale). See [NEWS.md](NEWS.md) for the details.

## Documentation

```r
help(package = "funresMech")
?funresMech
?run_app
```

## Citation

If you use funresMech in your research, please cite:

```bibtex
@article{NunezCampero2026,
  author = {Segundo Núñez-Campero},
  title = {funresMech: Mechanistic Functional Response Analysis using the Okuyama Model},
  year = {2026},
  note = {R package version 1.1.0},
  url = {https://github.com/Segon03/funresMech}
}

@article{Okuyama2012,
  author = {Okuyama, Toshinori},
  title = {A likelihood approach for functional response models},
  journal = {Biological Control},
  volume = {60},
  number = {2},
  pages = {103--107},
  year = {2012},
  doi = {10.1016/j.biocontrol.2011.10.008}
}

@article{Okuyama2026,
  author = {Okuyama, Toshinori},
  title = {Parametric Assumptions in Parasitoid Functional Response Analysis},
  journal = {Journal of Applied Entomology},
  year = {2026},
  doi = {10.1111/jen.70148}
}
```

## License

MIT License (see the `LICENSE` file). Copyright holder: Segundo Núñez-Campero.

## Contributing

Contributions are welcome. Please submit issues, feature requests or pull requests on GitHub.

## References

Okuyama, T. (2012). A likelihood approach for functional response models. *Biological Control*, 60(2), 103-107.

Okuyama, T. (2026). Parametric Assumptions in Parasitoid Functional Response Analysis. *Journal of Applied Entomology*.
