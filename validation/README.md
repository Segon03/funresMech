# funresMech: validation materials

Scripts, data, result tables, and figures of the verification and validation
of **funresMech 1.2.0**, reported in Núñez-Campero (2026), *funresMech: an R
package for mechanistic functional response analysis with simulation-based
likelihood*, and in its S1 Text. This folder is not part of the R package (it
is excluded from the package tarball through `.Rbuildignore`).

## What was evaluated

| Script | Content | S1 Text |
|---|---|---|
| `R/01_extract_forage_D2_D3.R` | Extraction of the published datasets D2 and D3 from the FoRAGE database | section 11 |
| `R/02_engine_tests.R` | Engine tests: anchor value of the supplementary code of the method, saturation stop, moments of the distributions | section 2 |
| `R/11_optimizer_comparison.R` | Four optimization strategies with equal budgets; profile against joint estimation of *z* | section 8 |
| `R/12_extended_validation.R` | Misspecified models, computational budget, and D2 and D3 with the final preset | sections 9 to 11 |
| `R/15_simulation_study.R` | Simulation study: 216 datasets, 288 fits, three designs, two presets | section 7 |
| `R/16_simulation_study_analysis.R` | Tests and additional analyses of the simulation study | section 7 |
| `R/17_parasitoid_vs_predator.R` | Parasitoid and predator processes in D2 and D3 | section 11 |
| `R/20_figures_and_tables.R` | Figures 1 to 3, S1 and S2 Figs, and the summary tables | all |
| `R/helpers.R` | Shared definitions: designs, scenarios, simulation of datasets, seeds | |

## Layout

| Folder | Content |
|---|---|
| `R/` | Numbered scripts. Run them from the root of the package repository, with funresMech installed. |
| `data/` | `D2.csv` and `D3.csv` (columns `dens`, host density, and `par`, hosts parasitized; *T* = 24 h), `forage_D2_D3_raw.csv` (with metadata), and the seed tables `seeds_datasets.csv` and `seeds_fits.csv`. |
| `tables/` | Result tables (CSV) of the runs reported in the article. |
| `figures/` | Figures of the article and of S1 Text (PNG; the TIFF files are written by `R/20_figures_and_tables.R`). |
| `results/` | Created by the scripts: one file per fit, so that an interrupted run resumes where it stopped. Not distributed. |

## Reproducing the results

```r
# From the root of the package repository:
install.packages(c("funresMech", "cmaes", "nloptr", "ggplot2"))
source("validation/R/02_engine_tests.R")                          # a few minutes
QUICK <- TRUE; source("validation/R/15_simulation_study.R")       # quick check of the installation
rm(list = ls()); source("validation/R/15_simulation_study.R")     # complete run (about one day on 23 cores)
source("validation/R/16_simulation_study_analysis.R")
source("validation/R/20_figures_and_tables.R")
```

Every script has a quick mode (`QUICK <- TRUE`) with very small budgets, to
check the installation in minutes before the complete runs. The complete runs
of scripts 11, 12, 15, and 17 took about four days on a workstation with 24
processor cores. Times per fit are reported in the result tables (`seconds`).

**Seeds.** Every simulated dataset and every fit uses a fixed seed, read from
`data/seeds_datasets.csv` (one base seed per scenario; the seed of dataset *j*
is base × 1000 + *j*) and `data/seeds_fits.csv` (one seed per fit). These
tables reproduce the seeds of the runs reported in the article. Scripts 11
and 02 use seeds written in the code.

**Code version.** The misspecification and budget fits and the fits of D2 and
D3 with the final preset (script 12) were run before the local refinement of
the profile of *z* and the shape-preserving interpolation of the interval
limits were added to the package. Their interval metrics can therefore differ
slightly in a new run; the intervals of D2 and D3 in
`tables/D2_D3_final_preset_fits.csv` were recomputed from the saved profiles
with the final interpolation (the linear limits are kept in the columns
`z_low_linear` and `z_high_linear`).

## Not included

The FoRAGE database (`FoRAGE_db_V5_Dec_20_2024_original_curves.csv`) is not
redistributed. Obtain it from Uiterwaal et al. (2022, *Ecology* 103: e3706)
and place it in `validation/data_raw/` to run `R/01_extract_forage_D2_D3.R`.
The extracted datasets are provided in `data/`.

## Main results

* With identical seeds, the engine reproduces the negative log-likelihood of
  the supplementary code of the method (132.04994217 for D3 at *z* = 2.5).
* In 288 fits to simulated data, the estimates of *z* and *h* were unbiased,
  *k* was overestimated, and there was no optimization failure. The 95%
  interval of *z* covered the true value in 91.0% of the fits (95% confidence
  interval 87.1% to 93.8%), and in 94.4% with 13 host densities.
* Differential evolution with common random numbers and a Nelder–Mead polish
  reached the best optimum within 0.5 NLL units in 80% of the runs with the
  exploratory budget and in all runs with the standard preset.
* Misspecified distributions of search or handling times had small effects;
  fitting the wrong process (parasitoid or predator) biased *h* and *z*.
