# ============================================================================
# 20_figures_and_tables.R
#
# Figures and summary tables of the manuscript and of S1 Text, computed from
# the result tables in Tables/ (produced by the validation scripts).
#
# Usage, either from the root of the package repository (tables in
# validation/tables, figures written to validation/figures):
#   source("validation/R/20_figures_and_tables.R")
# or from the folder of the submission that contains Tables/ and Figures/:
#   Rscript Scripts/20_figures_and_tables.R
# Requires: ggplot2 (>= 3.4.0), funresMech (>= 1.2.0, only for the
# interpolation of the profile interval in S2 Fig).
# ============================================================================

suppressMessages({ library(ggplot2) })
if (dir.exists(file.path("validation", "tables"))) {
  TAB <- file.path("validation", "tables"); FIG <- file.path("validation", "figures")
} else {
  TAB <- "Tables"; FIG <- "Figures"
}
dir.create(FIG, showWarnings = FALSE)
rd <- function(x) read.csv(file.path(TAB, x), stringsAsFactors = FALSE)
wr <- function(x, f) write.csv(x, file.path(TAB, f), row.names = FALSE)
wilson <- function(x) {
  n <- length(x); p <- mean(x); z <- 1.96
  c(p = p, lower = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    upper = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), n = n)
}
save_fig <- function(g, name, width, height) {
  ggsave(file.path(FIG, paste0(name, ".tif")), g, device = "tiff", width = width, height = height, dpi = 300,
         compression = "lzw", bg = "white")
  ggsave(file.path(FIG, paste0(name, ".png")), g, width = width, height = height, dpi = 300, bg = "white")
}
theme_set(theme_bw(base_size = 9) +
            theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
                  legend.position = "bottom"))
arm_labels <- c(core_exploratory = "8 or 5 densities, exploratory preset",
                core_standard    = "8 or 5 densities, standard preset",
                dense_design     = "13 densities, exploratory preset")
arm_cols <- c(core_exploratory = "#1F5F8B", core_standard = "#D9822B", dense_design = "#3A8F5C")

# ---- Fig 1: the simulated trial and the workflow -------------------------------
# (A) one simulated trial under each process; (B) workflow of an analysis.
set.seed(3)
trial <- function(process, x = 6, a = 0.08, z = 1, h = 1.2, T = 24) {
  t <- 0; out <- NULL; avail <- x; seen <- integer(0)
  while (t < T) {
    rate <- a * (if (process == "predator") avail else x)^z
    if (rate <= 0) { out <- rbind(out, data.frame(start = t, end = T, phase = "Search", host = NA)); break }
    s <- rgamma(1, shape = 1, rate = rate)
    out <- rbind(out, data.frame(start = t, end = min(T, t + s), phase = "Search", host = NA))
    t <- t + s; if (t >= T) break
    host <- if (process == "predator") setdiff(seq_len(x), seen)[sample.int(avail, 1)] else sample.int(x, 1)
    hd <- rlnorm(1, log(h) - 0.08, 0.4)
    out <- rbind(out, data.frame(start = t, end = min(T, t + hd), phase = "Handling", host = host))
    seen <- c(seen, host); if (process == "predator") avail <- avail - 1
    t <- t + hd
  }
  out$process <- process; out
}
tr <- rbind(trial("parasitoid"), trial("predator"))
tr$process <- factor(tr$process, c("parasitoid", "predator"),
                     c("Parasitoid process: attacked hosts remain available",
                       "Predator process: attacked hosts are removed"))
lab <- do.call(rbind, lapply(split(tr, tr$process), function(q) {
  hh <- q[q$phase == "Handling", ]
  if (!nrow(hh)) return(NULL)
  hh$label <- ifelse(duplicated(hh$host), paste0(hh$host, "*"), hh$host)
  hh }))
nresp <- vapply(split(tr, tr$process), function(q) length(unique(na.omit(q$host))), 0)
gA <- ggplot(tr) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = phase), colour = "white", linewidth = 0.2) +
  geom_text(data = lab, aes(x = (start + end) / 2, y = 0.5, label = label), size = 2.4, colour = "white") +
  facet_wrap(~process, ncol = 1) +
  scale_fill_manual(values = c(Search = "#9CB7CC", Handling = "#1F5F8B"), name = NULL) +
  scale_x_continuous("Time within the trial (h)", limits = c(0, 24), expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL) +
  geom_text(data = data.frame(process = factor(levels(tr$process), levels(tr$process)),
                              txt = paste("Response =", nresp, "distinct hosts")),
            aes(x = 23.8, y = 1.25, label = txt), hjust = 1, size = 2.6) +
  coord_cartesian(ylim = c(0, 1.45), clip = "off") +
  labs(title = "A. One simulated trial (6 hosts, T = 24 h)",
       caption = "Numbers: host attacked; * host attacked again (parasitoid process only).") +
  theme(plot.title = element_text(size = 9, face = "bold"), plot.caption = element_text(size = 7))
steps <- data.frame(
  id = 1:6,
  title = c("1. Data", "2. Screening", "3. Model", "4. Fit", "5. Reliability", "6. Report"),
  text = c("one row per trial:\ndensity and count",
           "screen_outliers()\nflag, never remove",
           "funres_model()\nprocess, distributions,\nfixed parameters",
           "fit_funres()\npreset, profile of z,\n95% interval, AIC",
           "verify_fit()\nreliability()\nact on every warning",
           "estimates and\ninterval of z,\nchecks, seed,\nversion"))
steps$x <- steps$id
gB <- ggplot(steps) +
  geom_segment(data = steps[-6, ], aes(x = x + 0.42, xend = x + 0.58, y = 0.5, yend = 0.5),
               arrow = arrow(length = unit(1.5, "mm")), colour = "grey40") +
  geom_tile(aes(x = x, y = 0.5), width = 0.84, height = 0.9, fill = "grey96", colour = "grey55") +
  geom_text(aes(x = x, y = 0.82, label = title), size = 2.6, fontface = "bold") +
  geom_text(aes(x = x, y = 0.42, label = text), size = 2.2, lineheight = 0.9) +
  scale_x_continuous(NULL, breaks = NULL, limits = c(0.5, 6.5), expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, limits = c(0, 1)) +
  labs(title = "B. Workflow of an analysis (all steps are also available in run_app())") +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(size = 9, face = "bold"))
if (requireNamespace("patchwork", quietly = TRUE)) {
  g1 <- patchwork::wrap_plots(gA, gB, ncol = 1, heights = c(1.6, 1))
} else {
  g1 <- gridExtra::arrangeGrob(gA, gB, ncol = 1, heights = c(1.6, 1))
}
save_fig(g1, "Fig1", 7.0, 4.6)

# ---- Fig 2 and Table 1: simulation study ---------------------------------------
G <- rd("simulation_study_fits.csv")
G$excludes_1 <- (!is.na(G$z_low) & G$z_low > 1) | (!is.na(G$z_high) & G$z_high < 1)  # both directions
G$true_z <- factor(G$z)
cov <- do.call(rbind, lapply(split(G, list(G$arm, G$z), drop = TRUE), function(q) {
  w <- wilson(q$covers)
  data.frame(arm = q$arm[1], z = q$z[1], n = nrow(q), coverage = w[["p"]], lower = w[["lower"]],
             upper = w[["upper"]])
}))
cov$arm <- factor(cov$arm, names(arm_labels))
g2a <- ggplot(cov, aes(factor(z), coverage, colour = arm)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(0.5), size = 0.3) +
  scale_colour_manual(values = arm_cols, labels = arm_labels, name = NULL) +
  scale_y_continuous("Coverage of the 95% interval of z", limits = c(0.5, 1)) +
  labs(x = "True z", title = "B. Interval coverage") +
  guides(colour = guide_legend(ncol = 1)) +
  theme(plot.title = element_text(size = 9, face = "bold"))
G$k_lab <- factor(paste("k =", G$k), paste("k =", c(0.2, 1, 5)))
g2b <- ggplot(G, aes(true_z, z_hat)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6, fill = "grey95") +
  geom_point(data = unique(G[, c("true_z", "z", "k_lab")]), aes(true_z, z), shape = 4, colour = "red3",
             size = 2) +
  facet_wrap(~k_lab) +
  labs(x = "True z", y = "Estimated z", title = "A. Estimates of z (all 288 fits)") +
  theme(plot.title = element_text(size = 9, face = "bold"))
g2 <- if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(g2b, g2a, widths = c(1.4, 1)) else
  gridExtra::arrangeGrob(g2b, g2a, ncol = 2, widths = c(1.4, 1))
save_fig(g2, "Fig2", 7.0, 3.6)

t1 <- do.call(rbind, lapply(split(G, list(G$arm, G$z), drop = TRUE), function(q) {
  w <- wilson(q$covers); pw <- wilson(q$excludes_1)
  data.frame(arm = q$arm[1], true_z = q$z[1], fits = nrow(q),
             median_z_hat = median(q$z_hat), median_abs_error = median(abs(q$z_hat - q$z)),
             coverage = w[["p"]], coverage_lower = w[["lower"]], coverage_upper = w[["upper"]],
             interval_excludes_1 = pw[["p"]], median_width = median(q$z_high - q$z_low, na.rm = TRUE),
             open_intervals = mean(q$open), median_h_ratio = median(q$h_hat / (q$h_over_T * 24)),
             median_k_ratio = median(q$k_hat / q$k), optimization_failures = sum(q$optimization_failure),
             mean_hours_one_core = mean(q$seconds) / 3600)
}))
t1 <- t1[order(match(t1$arm, names(arm_labels)), t1$true_z), ]
wr(t1, "Table2_simulation_study_summary.csv")
overall <- do.call(rbind, lapply(c(split(G, G$arm), list(all = G)), function(q) {
  w <- wilson(q$covers)
  data.frame(fits = nrow(q), coverage = w[["p"]], lower = w[["lower"]], upper = w[["upper"]])
}))
overall$arm <- rownames(overall)
wr(overall, "simulation_study_overall_coverage.csv")
fac <- G[G$variant == "exploratory", ]
by_factor <- rbind(
  do.call(rbind, lapply(split(fac, fac$k), function(q) { w <- wilson(q$covers)
    data.frame(factor = "k", level = q$k[1], fits = nrow(q), coverage = w[["p"]], lower = w[["lower"]], upper = w[["upper"]]) })),
  do.call(rbind, lapply(split(fac, fac$design), function(q) { w <- wilson(q$covers)
    data.frame(factor = "design", level = q$design[1], fits = nrow(q), coverage = w[["p"]], lower = w[["lower"]], upper = w[["upper"]]) })))
wr(by_factor, "simulation_study_coverage_by_factor.csv")

# Paired comparison of the presets on the same datasets
e <- G[G$arm == "core_exploratory", ]; s <- G[G$arm == "core_standard", ]
m <- merge(s, e, by = c("scenario", "dataset"), suffixes = c("_standard", "_exploratory"))
wr(data.frame(paired_fits = nrow(m), coverage_exploratory = mean(m$covers_exploratory),
              coverage_standard = mean(m$covers_standard),
              discordant = sum(m$covers_exploratory != m$covers_standard),
              precise_nll_better_by_0.5_with_standard = mean(m$nll_precise_hat_standard < m$nll_precise_hat_exploratory - 0.5),
              hours_exploratory = mean(m$seconds_exploratory) / 3600, hours_standard = mean(m$seconds_standard) / 3600),
   "preset_paired_comparison.csv")

# ---- Table 2 data: calibration of the reliability checks ------------------------
E <- rd("simulation_study_checks.csv")
key <- paste(G$scenario, G$dataset, G$arm)
E$good <- E$covers & !E$optimization_failure & abs(E$z_error) <= 1
cal <- do.call(rbind, lapply(split(E, list(E$variant, E$check), drop = TRUE), function(q)
  data.frame(preset = q$variant[1], check = q$check[1], fits = nrow(q),
             alerts = mean(q$status != "ok"), failures = mean(q$status == "fail"),
             alerts_in_good_fits = mean(q$status[q$good] != "ok"))))
wr(cal[order(cal$preset, cal$check), ], "reliability_checks_calibration.csv")

# ---- Fig 3: misspecified models -------------------------------------------------
B <- rd("misspecification_fits.csv")
B$z_error <- B$z_hat - B$z
gm <- c(correct = "Correct model", search_lognormal = "Lognormal search", search_weibull = "Weibull search",
        handling_gamma = "Gamma handling", predator = "Predator process",
        heterogeneity_a = "Individual variation in a")
B$gen_lab <- factor(gm[B$generating_model], gm)
B$fit_lab <- factor(ifelse(B$variant == "default_fit", "Default model fitted", "Generating model fitted"),
                    c("Default model fitted", "Generating model fitted"))
B$z_lab <- factor(paste("True z =", B$z))
g3 <- ggplot(B, aes(gen_lab, z_error, colour = fit_lab)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), size = 0.9, alpha = 0.7) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, position = position_dodge(0.6)) +
  facet_wrap(~z_lab, ncol = 1) +
  coord_cartesian(ylim = c(-1.5, 3.5)) +
  scale_colour_manual(values = c("#7A7A7A", "#1F5F8B"), name = NULL) +
  labs(x = "Process that generated the data", y = "Estimated z minus true z",
       caption = paste0("Bars: medians. Not shown (outside the vertical range): ",
                        sum(B$z_error > 3.5 | B$z_error < -1.5), " fit(s).")) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_fig(g3, "Fig3", 7.0, 4.6)
mis <- aggregate(cbind(z_error, covers, excludes_1, h_ratio = h_hat / (h_over_T * 24)) ~ generating_model + variant + z,
                 B, function(v) median(as.numeric(v)))
names(mis)[names(mis) == "covers"] <- "median_covers"
wr(do.call(rbind, lapply(split(B, list(B$generating_model, B$variant, B$z), drop = TRUE), function(q)
  data.frame(generating_model = q$generating_model[1], fit = q$variant[1], true_z = q$z[1], fits = nrow(q),
             median_z_error = median(q$z_error), coverage = mean(q$covers), interval_excludes_1 = mean(q$excludes_1),
             median_h_ratio = median(q$h_hat / (q$h_over_T * 24))))),
   "misspecification_summary.csv")
w <- merge(B[B$variant == "default_fit", c("scenario", "dataset", "nll")],
           B[B$variant == "correct_fit", c("scenario", "dataset", "nll", "generating_model")],
           by = c("scenario", "dataset"), suffixes = c("_default", "_correct"))
w$dAIC <- 2 * (w$nll_default - w$nll_correct)
wr(do.call(rbind, lapply(split(w, w$generating_model), function(q)
  data.frame(generating_model = q$generating_model[1], fits = nrow(q), median_dAIC = median(q$dAIC),
             correct_model_better = mean(q$dAIC > 0), correct_model_better_by_2 = mean(q$dAIC > 2)))),
   "misspecification_model_selection.csv")

# ---- S1 Fig: optimizer comparison -----------------------------------------------
O <- rd("optimizer_comparison_fits.csv")
st <- c(DEoptim = "Differential evolution", CMAES = "CMA-ES", CRS2 = "CRS2",
        DEoptim_CRN_polish = "Differential evolution,\ncommon random numbers\nand Nelder-Mead polish")
O$strategy_lab <- factor(ifelse(O$strategy %in% names(st), st[O$strategy], O$strategy), unique(c(st, O$strategy)))
gS1 <- ggplot(O, aes(strategy_lab, pmax(gap, 0.01))) +
  geom_boxplot(outlier.size = 0.6, fill = "grey95") +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  scale_y_log10("Gap to the best optimum found (NLL units)") +
  labs(x = NULL) + coord_flip()
save_fig(gS1, "S1_Fig", 6.0, 3.2)
wr(do.call(rbind, lapply(split(O, O$strategy), function(q)
  data.frame(strategy = q$strategy[1], fits = nrow(q), success_0.5 = mean(q$gap < 0.5), success_1 = mean(q$gap < 1),
             success_2 = mean(q$gap < 2), mean_gap = mean(q$gap), p95_gap = unname(quantile(q$gap, 0.95)),
             max_gap = max(q$gap), on_bound = mean(q$on_bound), seconds = median(q$seconds)))),
   "optimizer_comparison_summary.csv")
J <- rd("joint_vs_profile_fits.csv")
wr(do.call(rbind, lapply(split(J, J$strategy), function(q)
  data.frame(strategy = q$strategy[1], fits = nrow(q), coverage = mean(q$covers), open = mean(q$open),
             median_width = median(q$z_high - q$z_low, na.rm = TRUE),
             excludes_1_when_z_1 = mean(q$excludes_1[q$z_true == 1]),
             excludes_1_when_z_above_1 = mean(q$excludes_1[q$z_true > 1]),
             core_hours = sum(q$seconds) / 3600))),
   "joint_vs_profile_summary.csv")
Cb <- rd("budget_fits.csv")
wr(do.call(rbind, lapply(split(Cb, Cb$variant), function(q)
  data.frame(variant = q$variant[1], fits = nrow(q), success_0.5 = mean(q$success_0.5), success_1 = mean(q$success_1),
             mean_gap = mean(q$gap), median_hours_one_core = median(q$seconds) / 3600))),
   "budget_summary.csv")

# ---- S2 Fig: D2 and D3 ----------------------------------------------------------
FP <- rd("D2_D3_final_preset_profiles.csv")
PP <- rd("process_comparison_profiles.csv")
dl <- c(D2 = "D2", D2_without_atypical = "D2 without the atypical trial", D3 = "D3")
a1 <- do.call(rbind, lapply(split(FP, FP$data), function(q) { q$dnll <- q$nll - min(q$nll); q }))
a1$panel <- "Final preset, parasitoid process"; a1$series <- "parasitoid, final preset"
a2 <- do.call(rbind, lapply(split(PP, list(PP$data, PP$process, PP$seed_run), drop = TRUE), function(q) {
  q$dnll <- q$nll - min(PP$nll[PP$data == q$data[1] & PP$seed_run == q$seed_run[1]]); q }))
a2$panel <- "Exploratory preset, parasitoid and predator processes"
a2$series <- paste(a2$process, "run", a2$seed_run)
cols2 <- c("parasitoid, final preset" = "#0B2E4A", "parasitoid run 1" = "#1F5F8B", "parasitoid run 2" = "#6FA3C7",
           "predator run 1" = "#C0392B", "predator run 2" = "#E59866")
A2 <- rbind(a1[, c("data", "z", "dnll", "panel", "series")], a2[, c("data", "z", "dnll", "panel", "series")])
A2$data <- factor(dl[A2$data], dl)
gS2 <- ggplot(A2[A2$z <= 10, ], aes(z, dnll, colour = series)) +
  geom_hline(yintercept = 1.92, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 0.4) + geom_point(size = 0.7) +
  facet_grid(panel ~ data, labeller = label_wrap_gen(24)) +
  coord_cartesian(ylim = c(0, 12)) +
  scale_colour_manual(values = cols2, name = NULL) +
  labs(x = "Density-scaling exponent z", y = "NLL minus the minimum")
save_fig(gS2, "S2_Fig", 7.0, 4.6)
cat("Figures and tables written.\n")
