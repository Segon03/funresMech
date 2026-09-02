# funresMech: Mechanistic Functional Response Analysis using the Okuyama Model

[![CRAN status](https://www.r-pkg.org/badges/version/funresMech)](https://CRAN.R-project.org/package=funresMech)
[![R-CMD-check](https://github.com/Segon03/funresMech/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Segon03/funresMech/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**funresMech** implements the mechanistic, stochastic functional response model proposed by Okuyama (2012) and extended to parasitoids in Okuyama (2026). Unlike traditional approaches that rely on heuristic distributions (binomial or beta-binomial), this package simulates the underlying search-encounter-handling process to generate the probability distribution of the data, providing a more flexible and mechanistically sound framework for functional response analysis.

The package includes:
- An interactive **Shiny application** for data exploration and model fitting.
- **Maximum likelihood estimation** using a simulation-based likelihood.
- **Likelihood profiles** for the density-scaling exponent \\(z\\).
- **Model comparison** via AIC between full (\\(z\\) free) and restricted (\\(z = 1\\)) models.
- **Comprehensive diagnostic plots** including stochastic curves, histograms, density plots, boxplots, violins, and fan plots.

## Installation

### From GitHub (development version)
```r
# Install from GitHub using pak (recommended)
install.packages("pak")
pak::pkg_install("Segon03/funresMech")

# Or using devtools (legacy)
install.packages("devtools")
devtools::install_github("Segon03/funresMech")
From CRAN (stable version, once published)
r
install.packages("funresMech")
Basic Usage
Launch the Shiny App
r
library(funresMech)
run_app()
This opens the interactive application where you can:

Upload your dataset (CSV format).

Select columns for species, host density, and parasitism.

Configure advanced settings (simulation parameters, optimization options).

Run the analysis and explore results interactively.

Programmatic Usage (Advanced)
r
# Load the package
library(funresMech)

# Prepare your data (example format)
data <- data.frame(
  species = rep("Species_A", 30),
  dens = rep(c(10, 20, 40, 80, 160), each = 6),
  par = c(2, 3, 5, 8, 12, ...)  # Your data
)

# Fit the model (internal functions)
# See package documentation for details
Features
Mechanistic simulation: Search times follow a Gamma distribution; handling times follow a Lognormal distribution.

Stochastic likelihood: The probability distribution of parasitism is generated through repeated simulations.

Flexible density scaling: The exponent \(z\) allows emergence of Type I, II, III-like responses.

Uncertainty quantification: Confidence intervals for \(z\) via profile likelihood.

Interactive visualization: Dynamic plots with plotly for exploring results.

Comprehensive reporting: Generate HTML reports summarizing all analyses.

Documentation
Full documentation is available within the package:

r
# View package documentation
help(package = "funresMech")

# Get help for specific functions
?run_app
Citation
If you use funresMech in your research, please cite:

bibtex
@article{NunezCampero2026,
  author = {Segundo Núñez-Campero},
  title = {funresMech: Mechanistic Functional Response Analysis using the Okuyama Model},
  year = {2026},
  note = {R package version 1.0.0},
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
License
This package is distributed under the MIT License:

YEAR: 2026

COPYRIGHT HOLDER: Segundo Núñez-Campero

For more details, see the LICENSE file.

Contributing
Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests on GitHub.

References
Okuyama, T. (2012). A likelihood approach for functional response models. Biological Control, 60(2), 103–107.

Okuyama, T. (2026). Parametric Assumptions in Parasitoid Functional Response Analysis. Journal of Applied Entomology.