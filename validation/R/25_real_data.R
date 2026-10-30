# ============================================================================
# 25_real_data.R
#
# Published functional-response data analysed with the simulation-based
# method and with the conventional methods, following the selection protocol
# fixed before any fit (project document PROTOCOLO_seleccion_datos_reales,
# 2026-10-30).
#
# Datasets
#   FoRAGE V5 (DeLong and Uiterwaal, KNB, doi:10.5063/F1RX99KB): every curve
#   that meets all inclusion criteria
#     I1 one value per consumer and trial ("Rate is" = Point)
#     I2 taken from tables, supplements or shared data (not digitized from
#        figures, which loses overlapping points)
#     I3 integer counts not larger than the density offered
#     I4 one consumer per arena and no interference treatment
#     I5 known process; prey not replaced (the package has no process with
#        replacement): predator, or parasitoid with hosts left in the arena
#     I6 known trial duration
#     I7 at least 5 densities and 30 trials
#   Other published data that meet the same criteria (validation/data/real):
#     D'Auro et al. (2020), Pseudapanteles dignus on Tuta absoluta (dissected hosts)
#     Nunez-Campero et al. (2025), four parasitoids of Anastrepha fraterculus
#   Excluded: frair's gammarus (treatments pooled) and bythotrephes (24 trials
#   per group); D2 and D3 (digitized from a figure, criterion I2).
#
# Fits of each dataset
#   mechanistic, Gamma handling       funres_model(handling = "gamma", process)
#   mechanistic, lognormal handling   funres_model(handling = "lognormal", process)
#   conventional: holling, rogers, eq2 x binomial, beta-binomial, least squares
# The mechanistic fits use the exploratory settings followed by verify_fit()
# (screening); the dataset chosen for the illustration is refitted later with
# the default settings.
#
# Usage, from the root of the package repository, with the FoRAGE files in
# validation/data_raw/ (see validation/data/real/README.md):
#   QUICK <- TRUE; source("validation/R/25_real_data.R")   # 3 datasets, tiny settings
#   rm(list = ls()); source("validation/R/25_real_data.R")
# Optional: FAMILIES <- "gamma" to fit only Gamma handling.
# The run can be interrupted and resumed. Results:
#   validation/tables/real_datasets.csv     inventory with the reason of each exclusion
#   validation/tables/real_fits.csv         one row per dataset and method
#   validation/tables/real_comparison.csv   one row per dataset: AIC and conclusions by method
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")
source(file.path("validation", "R", "helpers.R"))
source(file.path("validation", "R", "conventional.R"))

if (!exists("QUICK")) QUICK <- FALSE
if (!exists("FAMILIES")) FAMILIES <- c("gamma", "lognormal")
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RAW <- file.path("validation", "data_raw")
DIR_RES <- file.path("validation", "results", if (QUICK) "real_quick" else "real")
DIR_TAB <- file.path("validation", "tables")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE)
SUFFIX <- if (QUICK) "_quick" else ""
GRID <- if (QUICK) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.5)
N_PRECISE <- if (QUICK) 2000L else 20000L
control_for <- function(...) {
  if (QUICK) effort_control("exploratory", n_sim = 200, itermax = 5, NP = 12, n_reeval = 2, ...)
  else effort_control("exploratory", ...)
}

# ---- FoRAGE: screening with the inclusion criteria --------------------------------
f_curves <- file.path(DIR_RAW, "FoRAGE_db_V5_Dec_20_2024_original_curves.csv")
f_meta   <- file.path(DIR_RAW, "FoRAGE_db_V5_Dec_20_2024_sources_and_meta.csv")
if (!file.exists(f_curves) || !file.exists(f_meta))
  stop("Place the FoRAGE V5 files in ", DIR_RAW, " (see validation/data/real/README.md).")
cu <- utils::read.csv(f_curves, check.names = FALSE, stringsAsFactors = FALSE, encoding = "latin1")
names(cu)[1] <- "id"
me <- utils::read.csv(f_meta, check.names = FALSE, stringsAsFactors = FALSE, encoding = "latin1")
me <- me[, !duplicated(names(me))]; names(me)[1] <- "id"
me <- me[!is.na(me$id), ]

screen_curve <- function(g) {
  x <- suppressWarnings(as.numeric(g[["Original x"]])); y <- suppressWarnings(as.numeric(g[["Original y"]]))
  m <- me[me$id == g$id[1], ]
  Tn <- suppressWarnings(as.numeric(g[["Trial duration (h)"]][1]))
  int <- all(is.finite(x) & is.finite(y)) && all(abs(y - round(y)) < 1e-6) && all(abs(x - round(x)) < 1e-6)
  reason <- c(
    I1 = g[["Rate is"]][1] != "Point",
    I2 = grepl("^fig", g[["Data origin"]][1], ignore.case = TRUE),
    I3 = !int || any(round(y) > round(x)),
    I4 = !nrow(m) || !(is.na(m[["Pred per arena"]]) || m[["Pred per arena"]] %in% c(1, "1")) ||
      identical(m[["Interference"]], "Y"),
    I5 = !nrow(m) || !identical(m[["Prey replaced?"]], "N"),
    I6 = !is.finite(Tn) || Tn <= 0,
    I7 = length(unique(round(x))) < 5 || length(x) < 30)
  data.frame(source = "FoRAGE", dataset = paste0("FoRAGE_", g$id[1]), study = if (nrow(m)) m$Source else NA,
             consumer = if (nrow(m)) m[["Predator scientific name"]] else g$Predator[1],
             resource = if (nrow(m)) m[["Prey scientific name"]] else g$Prey[1],
             process = if (nrow(m) && identical(m[["Predation or Parasitism"]], "Parasitoid")) "parasitoid" else "predator",
             T = Tn, n_trials = nrow(g), n_densities = length(unique(round(x))),
             included = !any(reason), excluded_by = paste(names(reason)[reason], collapse = ","),
             stringsAsFactors = FALSE)
}
INV <- do.call(rbind, lapply(split(cu, cu$id), screen_curve))
data_of <- list()
for (i in which(INV$included)) {
  g <- cu[paste0("FoRAGE_", cu$id) == INV$dataset[i], ]
  data_of[[INV$dataset[i]]] <- data.frame(dens = round(as.numeric(g[["Original x"]])),
                                          par = round(as.numeric(g[["Original y"]])))
}

# ---- other published data ---------------------------------------------------------
dr <- utils::read.csv(file.path("validation", "data", "real", "DAuro2020_trials.csv"))
data_of[["DAuro2020"]] <- dr[, c("dens", "par")]
nc <- utils::read.csv(file.path("validation", "data", "real", "NunezCampero2025.csv"), stringsAsFactors = FALSE)
spp_nc <- c(Ch = "Coptera haywardi", Dc = "Doryctobracon crawfordi", Gp = "Ganaspis pelleranoi", Ob = "Opius bellus")
for (s in names(spp_nc)) data_of[[paste0("NunezCampero2025_", s)]] <- nc[nc$spp == s, c("dens", "par")]
other <- data.frame(source = "other", dataset = c("DAuro2020", paste0("NunezCampero2025_", names(spp_nc))),
                    study = c("D'Auro et al. 2020 (Zenodo 4105529)", rep("Nunez-Campero et al. 2025 Insects", 4)),
                    consumer = c("Pseudapanteles dignus", unname(spp_nc)),
                    resource = c("Tuta absoluta", rep("Anastrepha fraterculus", 4)),
                    process = "parasitoid", T = 24, stringsAsFactors = FALSE)
other$n_trials <- vapply(other$dataset, function(k) nrow(data_of[[k]]), 0)
other$n_densities <- vapply(other$dataset, function(k) length(unique(data_of[[k]]$dens)), 0)
other$included <- TRUE; other$excluded_by <- ""
INV <- rbind(INV, other)
utils::write.csv(INV, file.path(DIR_TAB, paste0("real_datasets", SUFFIX, ".csv")), row.names = FALSE)
SEL <- INV[INV$included, ]
cat(sprintf("Datasets included: %d (%d studies); excluded: %d.\n", nrow(SEL), length(unique(SEL$study)),
            sum(!INV$included)))
if (QUICK) SEL <- SEL[c(1, which(SEL$process == "parasitoid")[1], nrow(SEL)), ]

# ---- jobs -------------------------------------------------------------------------
jobs <- expand.grid(dataset = SEL$dataset, method = c(paste0("mech_", FAMILIES), "conventional"),
                    stringsAsFactors = FALSE)
jobs$file <- file.path(DIR_RES, sprintf("%s_%s.rds", jobs$dataset, jobs$method))
pending <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Fits: %d in total, %d pending; %d workers.\n", nrow(jobs), nrow(pending), WORKERS))

run_one <- function(job) {
  library(funresMech)
  future::plan(future::sequential)
  t0 <- proc.time()[["elapsed"]]
  info <- SEL[SEL$dataset == job$dataset, ]
  d <- data_of[[job$dataset]]
  set.seed(sum(utf8ToInt(paste(job$dataset, job$method))))
  if (job$method == "conventional") {
    res <- do.call(rbind, lapply(names(MEANS), function(m) do.call(rbind, lapply(c("binomial", "bb", "ls"),
      function(err) analyse(d, m, err, T = info$T)))))
    res$excludes_1 <- excludes_one(res$z_low, res$z_high)
    out <- list(job = job[, c("dataset", "method")], res = res, seconds = proc.time()[["elapsed"]] - t0)
  } else {
    fam <- sub("^mech_", "", job$method)
    m <- funres_model(handling = fam, process = info$process)
    f <- suppressWarnings(fit_funres(d, info$T, m, z_grid = GRID, time_unit = "h", control = control_for(bb = TRUE)))
    f <- verify_fit(f)
    res <- fit_row(f)
    # precise likelihood at the estimate (more simulations), with the T of this dataset
    res$nll_precise <- mean(vapply(1:5, function(i) funresMech:::.nll_sim(f$par, m, d, info$T, N_PRECISE), 0))
    out <- list(job = job[, c("dataset", "method")], res = res, checks = checks_long(f), fit = lighten(f),
                seconds = proc.time()[["elapsed"]] - t0)
  }
  saveRDS(out, job$file)
  TRUE
}

run_jobs(pending, run_one,
         globals = mget(c("run_one", "pending", "SEL", "data_of", "GRID", "QUICK", "control_for", "EFFORT",
                          "effort_control", "N_PRECISE", "lighten", "fit_row", "excludes_one", "checks_long",
                          "MEANS", "LO", "UP", "objective", "par_index", "fit_one", "profile_value", "analyse",
                          "lambertW_log", "mean_rogers", "mean_eq2", "mean_holling", "T_EXP")),
         workers = WORKERS, log_file = file.path(DIR_RES, "errors.log"))

# ---- tables -----------------------------------------------------------------------
O <- lapply(jobs$file[file.exists(jobs$file)], readRDS)
if (length(O)) {
  type_of <- function(excl, z) ifelse(excl, ifelse(z > 1, "z > 1", "z < 1"), "z = 1 not rejected")
  R <- do.call(rbind, lapply(O, function(o) {
    r <- o$res
    if (o$job$method == "conventional") {
      data.frame(dataset = o$job$dataset, method = paste(r$mean, r$error, sep = "_"), z_hat = r$z_hat,
                 z_low = r$z_low, z_high = r$z_high, a_hat = r$a_hat, h_hat = r$h_hat, k_hat = NA_real_,
                 rho_hat = r$rho_hat, nll = r$nll, n_par = r$n_par, reliability = NA_character_,
                 excludes_1 = r$excludes_1, seconds = o$seconds / nrow(r), stringsAsFactors = FALSE)
    } else {
      data.frame(dataset = o$job$dataset, method = o$job$method, z_hat = r$z_hat, z_low = r$z_low,
                 z_high = r$z_high, a_hat = r$a_hat, h_hat = r$h_hat, k_hat = r$k_hat, rho_hat = NA_real_,
                 nll = r$nll_precise, n_par = 5L, reliability = r$reliability, excludes_1 = r$excludes_1,
                 seconds = o$seconds, stringsAsFactors = FALSE)
    }
  }))
  R <- merge(SEL[, c("dataset", "source", "study", "consumer", "resource", "process", "T", "n_trials", "n_densities")], R)
  R$AIC <- ifelse(is.na(R$nll), NA_real_, 2 * R$nll + 2 * R$n_par)
  R$type <- type_of(R$excludes_1, R$z_hat)
  utils::write.csv(R, file.path(DIR_TAB, paste0("real_fits", SUFFIX, ".csv")), row.names = FALSE)

  # one row per dataset: the mechanistic fit (Gamma handling) against the conventional fits with the
  # mean function of the process (rogers for a predator, eq2 for a parasitoid)
  C <- do.call(rbind, lapply(split(R, R$dataset), function(g) {
    mean_ok <- if (g$process[1] == "predator") "rogers" else "eq2"
    pick <- function(m) g[g$method == m, ][1, ]
    mech <- pick(paste0("mech_", FAMILIES[1]))
    bb <- pick(paste0(mean_ok, "_bb")); bi <- pick(paste0(mean_ok, "_binomial")); ls <- pick(paste0(mean_ok, "_ls"))
    fr <- pick("rogers_bb")                                  # what frair users would fit
    E <- mech$a_hat * exp(mean(log(unique(data_of[[g$dataset[1]]]$dens))))^mech$z_hat * g$T[1]
    data.frame(g[1, c("dataset", "source", "study", "consumer", "resource", "process", "T", "n_trials", "n_densities")],
               mech_reliability = mech$reliability, z_mech = mech$z_hat, z_mech_low = mech$z_low,
               z_mech_high = mech$z_high, type_mech = mech$type, k_mech = mech$k_hat, h_over_T = mech$h_hat / g$T[1],
               E = E, rho_bb = bb$rho_hat,
               dAIC_bb = bb$AIC - mech$AIC, dAIC_binomial = bi$AIC - mech$AIC,
               dAIC_gamma_vs_lognormal = if (all(c("mech_gamma", "mech_lognormal") %in% g$method))
                 pick("mech_lognormal")$AIC - pick("mech_gamma")$AIC else NA_real_,
               type_bb = bb$type, type_binomial = bi$type, type_ls = ls$type, type_frair_rogers_bb = fr$type,
               differs_from_bb = mech$type != bb$type, differs_from_ls = mech$type != ls$type,
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  utils::write.csv(C, file.path(DIR_TAB, paste0("real_comparison", SUFFIX, ".csv")), row.names = FALSE)
  cat(sprintf("\nDatasets fitted: %d. Mechanistic better than the beta-binomial by > 2 AIC: %d; worse by > 2: %d.\n",
              nrow(C), sum(C$dAIC_bb > 2, na.rm = TRUE), sum(C$dAIC_bb < -2, na.rm = TRUE)))
  cat(sprintf("Type of response different from the beta-binomial: %d datasets.\n", sum(C$differs_from_bb, na.rm = TRUE)))
}
