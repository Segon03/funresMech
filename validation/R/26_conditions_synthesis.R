# ============================================================================
# 26_conditions_synthesis.R
#
# Under which conditions do the conventional methods fail, and under which do
# they perform like the simulation-based method? Synthesis of the 960
# simulated datasets fitted by 22_conventional_methods.R (parasitoid: 576
# datasets, z = 0.5 to 2.75; predator: 384 datasets, z = 0.5 to 2.75; k =
# 0.2, 1, 5; base scenario and the three extended-range blocks).
#
# Methods compared, for each process:
#   mechanistic         the simulation-based fit with the generating process
#   correct curve + BB  the mean function of the process (eq2 for a parasitoid,
#                       rogers for a predator) with beta-binomial errors
#   correct curve + binomial, correct curve + least squares
#   frair default       rogers with beta-binomial errors, the usual choice in
#                       published analyses (for a predator it is the correct curve)
#
# Outputs
#   validation/tables/conditions_cells.csv      metrics by process, method, k and z
#   validation/tables/conditions_blocks.csv     metrics by process, method, block and k
#   validation/tables/conditions_glm.csv        logistic models of coverage by method
#   validation/figures/conditions_map.png       coverage, correct conclusion and bias of z
#
# Usage, from the root of the package repository (seconds):
#   source("validation/R/26_conditions_synthesis.R")
# ============================================================================

DIR_TAB <- file.path("validation", "tables")
DIR_FIG <- file.path("validation", "figures")
dir.create(DIR_FIG, showWarnings = FALSE)

A <- utils::read.csv(file.path(DIR_TAB, "conventional_fits.csv"), stringsAsFactors = FALSE)
# the interval of z excludes 1 in either direction (also for files written before 2026-10-30)
A$excludes_1 <- (!is.na(A$z_low) & A$z_low > 1) | (!is.na(A$z_high) & A$z_high < 1)
A$type_hat <- ifelse(A$excludes_1, ifelse(A$z_hat > 1, "above 1", "below 1"), "1")
A$type_true <- ifelse(A$z > 1, "above 1", ifelse(A$z < 1, "below 1", "1"))
A$correct_type <- A$type_hat == A$type_true
A$block <- ifelse(A$s_over_h != 0.5, "s/h = 1.5", ifelse(A$h_over_T != 0.025, "h/T = 0.1",
                  ifelse(A$E != 20, "E = 5", "base")))

label_of <- function(method, process) {
  ok <- if (process == "predator") "rogers" else "eq2"
  switch(method,
         mechanistic = "Mechanistic",
         rogers_bb = if (process == "predator") "Correct curve + BB" else "frair default (Rogers + BB)",
         stats::setNames(c("Correct curve + BB", "Correct curve + binomial", "Correct curve + LS"),
                         paste0(ok, c("_bb", "_binomial", "_ls")))[method])
}
A$label <- mapply(function(m, p) { l <- label_of(m, p); if (is.null(l) || is.na(l)) NA_character_ else unname(l) },
                  A$method, A$process)
A <- A[!is.na(A$label), ]
LEVELS <- c("Mechanistic", "Correct curve + BB", "Correct curve + binomial", "Correct curve + LS",
            "frair default (Rogers + BB)")
A$label <- factor(A$label, levels = LEVELS)

metrics <- function(g) data.frame(n = nrow(g), coverage = mean(g$covers), correct_type = mean(g$correct_type),
                                  median_z_error = stats::median(g$z_hat - g$z),
                                  median_h_ratio = stats::median(g$h_hat / g$h_true),
                                  excludes_1 = mean(g$excludes_1))
by_cells <- function(keys) {
  sp <- split(A, A[keys], drop = TRUE)
  out <- do.call(rbind, lapply(sp, function(g) cbind(g[1, keys, drop = FALSE], metrics(g))))
  rownames(out) <- NULL
  out[do.call(order, out[keys]), ]
}
CELLS <- by_cells(c("process", "label", "k", "z"))
BLOCKS <- by_cells(c("process", "label", "block", "k"))
utils::write.csv(CELLS, file.path(DIR_TAB, "conditions_cells.csv"), row.names = FALSE)
utils::write.csv(BLOCKS, file.path(DIR_TAB, "conditions_blocks.csv"), row.names = FALSE)

# Logistic models of coverage: which conditions lower it, for each method and process
GLM <- do.call(rbind, lapply(split(A, list(A$process, A$label), drop = TRUE), function(g) {
  m <- stats::glm(covers ~ factor(k) + factor(z) + block, family = stats::binomial, data = g)
  cf <- summary(m)$coefficients
  data.frame(process = g$process[1], method = as.character(g$label[1]), term = rownames(cf),
             odds_ratio = exp(cf[, 1]), p_value = cf[, 4], row.names = NULL)
}))
utils::write.csv(GLM, file.path(DIR_TAB, "conditions_glm.csv"), row.names = FALSE)

# ---- figure ---------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  long <- rbind(
    data.frame(CELLS[, c("process", "label", "k", "z")], metric = "Coverage of the 95% interval of z",
               value = CELLS$coverage, text = sprintf("%.0f", 100 * CELLS$coverage)),
    data.frame(CELLS[, c("process", "label", "k", "z")], metric = "Correct conclusion about z = 1",
               value = CELLS$correct_type, text = sprintf("%.0f", 100 * CELLS$correct_type)))
  long$label <- factor(long$label, levels = LEVELS)
  long$metric <- factor(long$metric, c("Coverage of the 95% interval of z", "Correct conclusion about z = 1"))
  long <- long[!(long$process == "predator" & long$label == "frair default (Rogers + BB)"), ]
  long$process <- factor(long$process, c("parasitoid", "predator"), c("Parasitoid process", "Predator process"))
  long$k <- factor(long$k, c(0.2, 1, 5), c("0.2 (irregular)", "1 (exponential)", "5 (regular)"))
  long$z <- factor(long$z)
  p <- ggplot(long, aes(z, k, fill = value)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = text), size = 2.6) +
    scale_fill_gradientn(colours = c("#b2182b", "#f4a582", "#f7f7f7", "#92c5de", "#2166ac"),
                         values = c(0, 0.5, 0.8, 0.9, 1), limits = c(0, 1), name = "Proportion") +
    facet_grid(process + metric ~ label, labeller = label_wrap_gen(18)) +
    labs(x = "True z", y = "Search dispersion k") +
    theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), legend.position = "bottom")
  ggsave(file.path(DIR_FIG, "conditions_map.png"), p, width = 10, height = 9, dpi = 300, bg = "white")
}

# ---- console summary -------------------------------------------------------------------
fmt <- function(x) sprintf("%.0f%%", 100 * x)
for (pr in c("parasitoid", "predator")) {
  cat(sprintf("\n%s process: coverage (correct conclusion) by k, pooled over z and blocks\n", pr))
  t <- by_cells(c("process", "label", "k")); t <- t[t$process == pr, ]
  print(stats::reshape(data.frame(label = t$label, k = t$k, v = paste0(fmt(t$coverage), " (", fmt(t$correct_type), ")")),
                       idvar = "label", timevar = "k", direction = "wide"), row.names = FALSE)
}
