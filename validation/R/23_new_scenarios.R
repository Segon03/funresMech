# ============================================================================
# 23_new_scenarios.R
#
# Simulation studies requested in the review of the manuscript, run with the
# same design and settings as 15_simulation_study.R and 18_extended_range.R
# (lower-effort settings, n_sim = 2000 and itermax = 50, followed by
# verify_fit(); z grid from 0.5 to 5 in steps of 0.5 with refinement and
# extension; 8 datasets per scenario). Each dataset is fitted with the model
# that generated it.
#
# STUDY <- "parasitoid_z05"  Parasitoid process with z = 0.5 (a response that
#   decelerates faster than type II), k = 0.2, 1 and 5: base scenario with 8,
#   5 and 13 densities, and the three extended-range blocks (s/h = 1.5,
#   h/T = 0.1, E = 5) with 8 densities. 18 scenarios, 144 fits.
# STUDY <- "predator"        Predator process (consumed prey removed),
#   z = 0.5, 1, 1.75 and 2.75, k = 0.2, 1 and 5: base scenario and the three
#   extended-range blocks, with 8 densities. 48 scenarios, 384 fits.
#
# Usage, from the root of the package repository with funresMech installed:
#   STUDY <- "parasitoid_z05"; QUICK <- TRUE; source("validation/R/23_new_scenarios.R")  # check
#   rm(list = ls()); STUDY <- "parasitoid_z05"; source("validation/R/23_new_scenarios.R")
#   rm(list = ls()); STUDY <- "predator";       source("validation/R/23_new_scenarios.R")
# Approximate cost: 130 core-hours for parasitoid_z05 (about 6 h on 23
# cores) and 350 core-hours for predator (about 15 h on 23 cores). The runs
# can be interrupted and resumed. Results:
#   validation/tables/<STUDY>_fits.csv, <STUDY>_checks.csv, <STUDY>_summary.csv
#   validation/data/simulated_datasets_<STUDY>.csv (the datasets, used by
#   22_conventional_methods.R)
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")
source(file.path("validation", "R", "helpers.R"))

if (!exists("STUDY")) stop('Set STUDY <- "parasitoid_z05" or STUDY <- "predator" first.')
STUDY <- match.arg(STUDY, c("parasitoid_z05", "predator"))
if (!exists("QUICK")) QUICK <- FALSE
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (QUICK) paste0(STUDY, "_quick") else STUDY)
DIR_TAB <- file.path("validation", "tables")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE); dir.create(DIR_TAB, showWarnings = FALSE)
SUFFIX <- if (QUICK) "_quick" else ""
GRID <- if (QUICK) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.5)
N_PRECISE <- if (QUICK) 2000L else 20000L

control_for <- function(level, ...) {
  if (QUICK) effort_control(level, n_sim = 200, itermax = 5, NP = 12, n_reeval = 2, ...)
  else effort_control(level, ...)
}

# ---- scenarios -------------------------------------------------------------------
K <- c(0.2, 1, 5)
BLOCKS <- list("s/h = 1.5" = list(s_over_h = 1.5), "h/T = 0.1" = list(h_over_T = 0.1), "E = 5" = list(E = 5))
make <- function(z_values, designs, gen) {
  base <- do.call(rbind, apply(expand.grid(z = z_values, k = K, design = designs, stringsAsFactors = FALSE), 1,
    function(r) scenario(z = as.numeric(r[["z"]]), k = as.numeric(r[["k"]]), design = r[["design"]],
                         generating_model = gen, block = "base")))
  ext <- do.call(rbind, lapply(names(BLOCKS), function(b)
    do.call(rbind, apply(expand.grid(z = z_values, k = K), 1, function(r)
      do.call(scenario, c(list(z = r[["z"]], k = r[["k"]], generating_model = gen, block = b), BLOCKS[[b]]))))))
  rbind(base, ext)
}
SCEN <- if (STUDY == "parasitoid_z05") make(0.5, c("large", "small", "dense"), "correct") else
  make(c(0.5, 1, 1.75, 2.75), "large", "predator")
SCEN$id <- scenario_id(SCEN)
stopifnot(!anyDuplicated(SCEN$id))

# Seeds: added once to the seed tables, then always read from them
add_seeds <- function() {
  sd <- utils::read.csv(file.path("validation", "data", "seeds_datasets.csv"), stringsAsFactors = FALSE)
  sf <- utils::read.csv(file.path("validation", "data", "seeds_fits.csv"), stringsAsFactors = FALSE)
  new <- setdiff(SCEN$id, sd$scenario)
  if (!length(new)) return(invisible(FALSE))
  base <- vapply(new, function(i) sum(utf8ToInt(i)) + 20000, 0)
  while (anyDuplicated(base) || any(base %in% sd$seed_base)) {
    j <- which(duplicated(base) | base %in% sd$seed_base)[1]; base[j] <- base[j] + 1
  }
  sd <- rbind(sd, data.frame(scenario = new, seed_base = unname(base)))
  sf <- rbind(sf, do.call(rbind, lapply(new, function(i)
    data.frame(script = "23", arm = STUDY, scenario = i, dataset = 1:8, variant = "exploratory",
               seed = sum(utf8ToInt(paste(STUDY, i))) * 7 + 13 + (1:8) * 101))))
  utils::write.csv(sd, file.path("validation", "data", "seeds_datasets.csv"), row.names = FALSE)
  utils::write.csv(sf, file.path("validation", "data", "seeds_fits.csv"), row.names = FALSE)
  invisible(TRUE)
}
if (add_seeds()) {                                  # reload the tables with the new seeds
  SEEDS_DATA <- utils::read.csv(file.path("validation", "data", "seeds_datasets.csv"), stringsAsFactors = FALSE)
  SEEDS_FIT  <- utils::read.csv(file.path("validation", "data", "seeds_fits.csv"), stringsAsFactors = FALSE)
}
if (QUICK) SCEN <- SCEN[unique(round(seq(1, nrow(SCEN), length.out = 3))), ]

jobs <- do.call(rbind, lapply(SCEN$id, function(id)
  data.frame(arm = STUDY, scenario = id, dataset = seq_len(if (QUICK) 1 else 8), variant = "exploratory")))
jobs$file <- file.path(DIR_RES, sprintf("%s_d%d.rds", jobs$scenario, jobs$dataset))
pending <- jobs[!file.exists(jobs$file), ]
cat(sprintf("%s: %d fits in total, %d pending; %d workers.\n", STUDY, nrow(jobs), nrow(pending), WORKERS))

run_one <- function(job) {
  library(funresMech)
  future::plan(future::sequential)
  t0 <- proc.time()[["elapsed"]]
  e <- SCEN[SCEN$id == job$scenario, ]
  d <- simulate_dataset(e, job$dataset); p0 <- true_parameters(e); m <- generating_model(e)
  set.seed(fit_seed("23", job$arm, job$scenario, job$dataset, job$variant))
  f <- suppressWarnings(fit_funres(d, T_EXP, m, z_grid = GRID, control = control_for(job$variant, bb = TRUE)))
  f <- verify_fit(f)
  res <- fit_row(f, e$z)
  res$nll_precise_hat <- precise_nll(f$par, m, d)
  res$nll_precise_true <- precise_nll(p0, m, d)
  out <- list(job = job[, c("arm", "scenario", "dataset", "variant")], res = res, checks = checks_long(f),
              fit = lighten(f), seconds = proc.time()[["elapsed"]] - t0)
  saveRDS(out, job$file)
  TRUE
}

run_jobs(pending, run_one,
         globals = mget(c("run_one", "pending", "SCEN", "DESIGNS", "T_EXP", "GRID", "N_PRECISE", "QUICK",
                          "control_for", "EFFORT", "effort_control", "true_parameters", "simulate_dataset",
                          "generating_model", "data_seed", "fit_seed", "SEEDS_DATA", "SEEDS_FIT", "precise_nll",
                          "lighten", "fit_row", "checks_long")),
         workers = WORKERS, log_file = file.path(DIR_RES, "errors.log"))

# ---- tables ------------------------------------------------------------------
O <- lapply(jobs$file[file.exists(jobs$file)], readRDS)
if (length(O)) {
  info <- c("z", "k", "s_over_h", "h_over_T", "E", "design", "block", "generating_model")
  G <- do.call(rbind, lapply(O, function(o)
    data.frame(o$job, o$res, seconds = o$seconds, SCEN[SCEN$id == o$job$scenario, info], row.names = NULL)))
  G$process <- ifelse(G$generating_model == "predator", "predator", "parasitoid")
  G$excludes_1 <- excludes_one(G$z_low, G$z_high)
  G$deficit <- G$nll_precise_hat - G$nll_precise_true
  G$optimization_failure <- G$deficit > 1
  G$open <- G$low_open | G$high_open
  G$z_error <- G$z_hat - G$z
  G$h_ratio <- G$h_hat / (G$h_over_T * T_EXP)
  G$k_ratio <- G$k_hat / G$k
  utils::write.csv(G, file.path(DIR_TAB, paste0(STUDY, "_fits", SUFFIX, ".csv")), row.names = FALSE)

  C <- do.call(rbind, lapply(O, function(o) {
    k <- merge(o$job, o$checks)
    k$covers <- o$res$covers
    k$optimization_failure <- (o$res$nll_precise_hat - o$res$nll_precise_true) > 1
    k }))
  utils::write.csv(C, file.path(DIR_TAB, paste0(STUDY, "_checks", SUFFIX, ".csv")), row.names = FALSE)

  # the simulated datasets, for the conventional analyses (22_conventional_methods.R)
  DS <- do.call(rbind, lapply(O, function(o)
    data.frame(source = STUDY, scenario = o$job$scenario, dataset = o$job$dataset, o$fit$data[, c("dens", "par")])))
  utils::write.csv(DS, file.path("validation", "data", paste0("simulated_datasets_", STUDY, SUFFIX, ".csv")),
                   row.names = FALSE)

  summarise <- function(g, label) {
    w <- wilson(g$covers)
    data.frame(group = label, n = nrow(g), coverage = w[["p"]], coverage_low = w[["lower"]],
               coverage_high = w[["upper"]],
               type_I_error = if (any(g$z == 1)) mean(g$excludes_1[g$z == 1]) else NA,
               power = if (any(g$z > 1)) mean(g$excludes_1[g$z > 1]) else NA,
               detects_z_below_1 = if (any(g$z < 1)) mean(g$excludes_1[g$z < 1]) else NA,
               open_interval = mean(g$open), optimization_failure = mean(g$optimization_failure),
               median_z_error = stats::median(g$z_error),
               median_h_ratio = stats::median(g$h_ratio), median_k_ratio = stats::median(g$k_ratio),
               median_hours = stats::median(g$seconds) / 3600)
  }
  S <- rbind(summarise(G, "all"),
             do.call(rbind, lapply(split(G, G$block), function(g) summarise(g, g$block[1]))),
             do.call(rbind, lapply(split(G, G$z), function(g) summarise(g, paste("z =", g$z[1])))),
             do.call(rbind, lapply(split(G, G$k), function(g) summarise(g, paste("k =", g$k[1])))))
  utils::write.csv(S, file.path(DIR_TAB, paste0(STUDY, "_summary", SUFFIX, ".csv")), row.names = FALSE)
  for (i in seq_len(nrow(S)))
    cat(sprintf("Coverage, %s: %.1f%% [%.1f, %.1f] (n = %d)\n", S$group[i], 100 * S$coverage[i],
                100 * S$coverage_low[i], 100 * S$coverage_high[i], S$n[i]))
}
