# ============================================================
# R/05_sintesis_y_figuras.R
#
# Combina los resultados de Ensayo 1 (tablas/ensayo1_*.csv) y
# Ensayo 2 (tablas/ensayo2_D2_D3_resultado.csv) en las figuras
# que van al paper.
#
# IMPORTANTE: correr 03_ensayo1_resumen.R y 04_ensayo2_datos_reales.R
# antes de esto.
# ============================================================

library(ggplot2)
library(dplyr)

dir.create("figuras", showWarnings = FALSE)

# ---- Ensayo 1: sesgo y cobertura de z vs rho -----------------------------
res1     <- read.csv(file.path("tablas", "ensayo1_ajustes_individuales.csv"))
resumen1 <- read.csv(file.path("tablas", "ensayo1_resumen_por_rho.csv"))

fig1 <- ggplot(res1, aes(x = factor(rho), y = z_hat)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_boxplot(fill = "#8FB6D9") +
  geom_jitter(width = 0.1, alpha = 0.4) +
  labs(x = expression(rho ~ "(sobre-dispersion beta-binomial)"),
       y = expression(hat(z)),
       title = "Ensayo 1: recuperacion de z (verdadero = 1) segun sobre-dispersion") +
  theme_bw(base_size = 13)
ggsave(file.path("figuras", "fig_ensayo1_z_por_rho.png"), fig1, width = 7, height = 5, dpi = 300)

fig2 <- ggplot(resumen1, aes(x = rho, y = cobertura_z)) +
  geom_line() + geom_point(size = 3) +
  geom_hline(yintercept = 0.95, linetype = "dotted") +
  ylim(0, 1) +
  labs(x = expression(rho), y = "Cobertura del IC95% de z",
       title = "Ensayo 1: cobertura del intervalo de confianza de z") +
  theme_bw(base_size = 13)
ggsave(file.path("figuras", "fig_ensayo1_cobertura.png"), fig2, width = 7, height = 5, dpi = 300)

# ---- Ensayo 2: perfiles de verosimilitud y forest plot de z --------------
perfiles <- lapply(c("D2", "D3"), function(lab) {
  r <- readRDS(file.path("results", sprintf("ensayo2_mech_%s.rds", lab)))
  data.frame(dataset = lab, z = r$z_grid, nll = r$nll_profile)
})
perfiles <- do.call(rbind, perfiles)

fig3 <- ggplot(perfiles, aes(x = z, y = nll)) +
  geom_line(color = "#2E5496") + geom_point(size = 1.5) +
  facet_wrap(~dataset, scales = "free_y") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  labs(x = "z", y = "Log-verosimilitud negativa",
       title = "Ensayo 2: perfil de verosimilitud de z (modelo mecanistico)") +
  theme_bw(base_size = 13)
ggsave(file.path("figuras", "fig_ensayo2_perfiles.png"), fig3, width = 8, height = 4.5, dpi = 300)

tabla2 <- read.csv(file.path("tablas", "ensayo2_D2_D3_resultado.csv"))
fig4 <- ggplot(tabla2, aes(x = dataset, y = z_hat_mech)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbar(aes(ymin = ci_low_mech, ymax = ci_high_mech), width = 0.15, linewidth = 1) +
  geom_point(size = 3, color = "#2E5496") +
  labs(x = NULL, y = "z estimado (IC95%)",
       title = "Ensayo 2: estimacion de z, modelo mecanistico") +
  theme_bw(base_size = 13)
ggsave(file.path("figuras", "fig_ensayo2_forest_z.png"), fig4, width = 5, height = 5, dpi = 300)

cat("Figuras guardadas en figuras/.\n")
