# funresMech: validation materials

Scripts, data, results and figures supporting the validation of **funresMech 1.1.0**
(C++/Rcpp engine) reported in Núñez-Campero (2026), *funresMech: An R Package for
Mechanistic Functional Response Analysis using the Okuyama Model*. This folder is **not**
part of the R package (it is excluded from the package tarball through `.Rbuildignore`).

Package: <https://CRAN.R-project.org/package=funresMech> (DOI 10.32614/CRAN.package.funresMech).

## What was validated

1. **Engine equivalence** (`R/00b_motor_rcpp_pruebas.R`): the C++ engine gives the same
   negative log-likelihood as the reference implementations for identical seeds.
2. **Parameter recovery** (Trial 1; `R/02c_ensayo1_rcpp.R`, `R/03b_ensayo1_rcpp_resumen.R`):
   3 scenarios x 20 simulated datasets (design of D3; T = 24 h).
3. **Reproduction of Okuyama (2026)** in two real datasets (Trial 2; `R/04h`-`R/04k`): D2 and D3,
   FoRAGE database v5 (*Apanteles myeloenta* on *Ectomyelois ceratoniae*;
   Kishani Farahani & Goldansaz 2013).
4. **Screening of atypical trials** (`R/06_diagnostico_atipicos.R`, `fun/diagnostico_datos.R`).
5. **Check of the released package** (`R/07_validacion_paquete.R`): the installed package
   (`funresMech:::negloglik_fixed_z`, `screen_outliers`, `fit_profile`) against the validated engine
   and against Okuyama (2026).

## Layout

| Folder | Content |
|---|---|
| `R/` | Numbered scripts (run from the project root, i.e. the folder with the `.Rproj`). `02`-`04g` are the archived pure-R-engine runs, kept as a record. |
| `fun/` | Project engine (`motor_okuyama.cpp`, `motor_rcpp.R`), beta-binomial model, profile utilities, outlier screening. |
| `data_clean/` | `D2.rds`, `D3.rds` and `forage_D2_D3_raw.csv` (host density `dens`, hosts parasitised `par`), extracted from FoRAGE v5. |
| `tablas/` | Result tables (CSV). Final tables: `ensayo1_rcpp_*`, `ensayo2_rcpp_*`, `diagnostico_*`, `paquete_*`. |
| `figuras/` | Figures. |
| `results/` | Fitted objects: `ensayo1_rcpp/` (simulated datasets and profile points), `ensayo2_rcpp/`, `paquete_1.1.0/` (`fit_D2.rds`, `fit_D3.rds`, `sessionInfo_*.txt`). |
| `docs/` | Audit (`AUDITORIA_validacion_2026-10-05.md`), engine notes (`MOTOR_RCPP_cambios_y_validacion.md`) and the project log (`LOG_proyecto_es.md`), in Spanish. |
| `report/` | `validacion_funresMech.Rmd` (reproducible report; not yet updated to the final tables). |

## Not included

The FoRAGE database file (`FoRAGE_db_V5_Dec_20_2024_original_curves.csv`, 16 MB) is not
redistributed; obtain it from Uiterwaal et al. (2022, *Ecology* 103(5): e3706) and place it in
`data_raw/` to re-run `R/01_extract_forage_d2_d3.R`. The extracted D2 and D3 are provided in `data_clean/`.
Temporary smoke-test runs (`*_PRUEBA`) and local backup files are also omitted.

## Reproducing the package check (`R/07_validacion_paquete.R`)

```r
remotes::install_github("Segon03/funresMech")   # or install.packages("funresMech")
# From the project root:
#   Rscript R/07_validacion_paquete.R prueba      # smoke test (minutes)
#   Rscript R/07_validacion_paquete.R D2 11       # D2 with 11 workers (about 14 h)
#   Rscript R/07_validacion_paquete.R D3 11       # D3 with 11 workers (about 14 h)
#   Rscript R/07_validacion_paquete.R comparar    # comparison tables only
```

Reference run (2026-10-06/07): funresMech 1.1.0 installed from GitHub `Segon03/funresMech` at commit
`2d12fe9`; R 4.6.1, Windows, Rcpp 1.1.2, DEoptim 2.2-8, future 1.75.0; 10,000 simulated trials per density,
population 40, 200 generations; fixed seeds (D2 20261006, D3 20261007).
The engine (`fun/`) requires a C++ toolchain (Rtools on Windows).
