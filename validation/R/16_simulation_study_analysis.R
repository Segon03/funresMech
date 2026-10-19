# ============================================================================
# 16_simulation_study_analysis.R
#
# Additional analyses of the simulation study (S1 Text, section 7). Reads the
# tables and fitted objects written by 15_simulation_study.R; refits nothing.
#
#   1. Checks of the run: every fit present, every interval contains its own
#      estimate, optimization failures, profile points with probability floor.
#   2. Coverage against the nominal 95% (exact binomial test), coverage by
#      design (Fisher's exact test), and the side of the misses.
#   3. Type I error and power by the interval and by DeltaAIC > 2.
#   4. Ability of each reliability check to predict a miss of the interval
#      (area under the ROC curve).
#   5. Post hoc: coverage of the profile interval for thresholds other than
#      1.92, computed on the same profiles.
#
# Usage: source("validation/R/16_simulation_study_analysis.R")
# Output: validation/tables/simulation_study_analysis_*.csv
# ============================================================================

suppressMessages(library(funresMech))
source(file.path("validation", "R", "helpers.R"))
DIR_TAB <- file.path("validation", "tables")
DIR_RES <- file.path("validation", "results", "simulation_study")
write_table <- function(x, name) write.csv(x, file.path(DIR_TAB, paste0("simulation_study_analysis_", name, ".csv")),
                                           row.names = FALSE)
auc <- function(x, y) {
  ok <- !is.na(x); x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(x); n1 <- sum(y); n0 <- sum(!y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

G <- read.csv(file.path(DIR_TAB, "simulation_study_fits.csv"), stringsAsFactors = FALSE)
C <- read.csv(file.path(DIR_TAB, "simulation_study_checks.csv"), stringsAsFactors = FALSE)
G$miss_side <- ifelse(G$covers, "covers", ifelse(!is.na(G$z_low) & G$z < G$z_low, "true z below the interval",
                                                 "true z above the interval"))

# ---- 1. Checks of the run ----------------------------------------------------
run_checks <- data.frame(
  check = c("fits", "interval contains its own estimate", "optimization failures",
            "fits with profile points on the probability floor"),
  value = c(nrow(G),
            with(G, all((is.na(z_low) | z_low <= z_hat + 1e-9) & (is.na(z_high) | z_hat <= z_high + 1e-9))),
            sum(G$optimization_failure), sum(G$floor_points > 0)))
write_table(run_checks, "run_checks"); print(run_checks)

# ---- 2. Coverage ---------------------------------------------------------------
bt <- stats::binom.test(sum(G$covers), nrow(G), p = 0.95)
expl <- G[G$variant == "exploratory", ]
fisher_design <- stats::fisher.test(table(expl$design, expl$covers))
coverage <- data.frame(
  measure = c("overall coverage", "exact binomial test against 0.95 (p)",
              "Fisher's exact test, coverage by design, exploratory fits (p)",
              "misses with the true z below the interval", "misses with the true z above the interval"),
  value = c(mean(G$covers), bt$p.value, fisher_design$p.value,
            sum(G$miss_side == "true z below the interval"), sum(G$miss_side == "true z above the interval")))
write_table(coverage, "coverage_tests"); print(coverage, digits = 3)

# ---- 3. Type I error, power, and bias ------------------------------------------
power <- do.call(rbind, lapply(split(G, list(G$arm, G$z), drop = TRUE), function(q) {
  w <- wilson(q$excludes_1)
  data.frame(arm = q$arm[1], true_z = q$z[1], fits = nrow(q),
             measure = if (q$z[1] == 1) "type I error" else "power",
             interval_excludes_1 = w[["p"]], lower = w[["lower"]], upper = w[["upper"]],
             dAIC_above_2 = mean(q$dAIC_z1 > 2), median_z_error = median(q$z_error),
             median_h_ratio = median(q$h_hat / (q$h_over_T * T_EXP)), median_k_ratio = median(q$k_hat / q$k),
             open_intervals = mean(q$open))
}))
write_table(power, "type_I_error_power_bias"); print(power, digits = 3)

# ---- 4. Reliability checks as predictors of a miss --------------------------------
C$good <- C$covers & !C$optimization_failure & abs(C$z_error) <= 1
pred <- do.call(rbind, lapply(split(C, list(C$variant, C$check), drop = TRUE), function(q)
  data.frame(preset = q$variant[1], check = q$check[1], fits = nrow(q),
             alerts = mean(q$status != "ok"), failures = mean(q$status == "fail"),
             alerts_in_good_fits = mean(q$status[q$good] != "ok"),
             auc_miss = auc(q$value, !q$covers))))
write_table(pred, "checks_prediction")
cat(sprintf("\nAUC of the checks for a miss: %.2f to %.2f\n", min(pred$auc_miss, na.rm = TRUE),
            max(pred$auc_miss, na.rm = TRUE)))

# ---- 5. Post hoc thresholds ------------------------------------------------------
interval_at <- function(prof, thr) {
  prof <- prof[prof$n_floor == 0 | is.na(prof$n_floor), ]
  prof <- prof[order(prof$z), ]
  f <- funresMech:::.pchip(prof$z, prof$nll)
  zz <- sort(unique(c(prof$z, seq(min(prof$z), max(prof$z), by = 0.01))))
  nn <- f(zz); i0 <- which.min(nn); cut <- min(nn) + thr
  lo <- if (any(nn[seq_len(i0)] >= cut)) { i <- max(which(nn[seq_len(i0)] >= cut))
    zz[i] + (zz[i + 1] - zz[i]) * (cut - nn[i]) / (nn[i + 1] - nn[i]) } else NA_real_
  hi <- if (any(nn[i0:length(nn)] >= cut)) { i <- i0 - 1 + min(which(nn[i0:length(nn)] >= cut))
    zz[i - 1] + (zz[i] - zz[i - 1]) * (cut - nn[i - 1]) / (nn[i] - nn[i - 1]) } else NA_real_
  c(lo, hi)
}
files <- list.files(DIR_RES, "\\.rds$", full.names = TRUE)
if (length(files)) {
  thresholds <- c(1.92, 2.25, 2.5, 2.75, 3.0, 3.5)
  cov_thr <- do.call(rbind, lapply(files, function(fl) {
    o <- readRDS(fl); z0 <- G$z[G$scenario == o$job$scenario][1]
    vapply(thresholds, function(t) {
      ci <- interval_at(o$fit$profile, t)
      (is.na(ci[1]) || ci[1] <= z0) && (is.na(ci[2]) || z0 <= ci[2])
    }, logical(1))
  }))
  post <- data.frame(threshold = thresholds, coverage = colMeans(cov_thr))
  write_table(post, "post_hoc_thresholds"); print(post, digits = 3)
}
