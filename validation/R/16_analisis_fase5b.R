# ============================================================================
# 16_analisis_fase5b.R
#
# Verificación e interpretación de la Fase 5b (15_confirmacion_fase5b.R).
# Lee validation/tablas/fase5b_*.csv, validation/results/fase5b/*.rds y, para
# las comparaciones pareadas, validation/results/fase5/A__*.rds.
# Escribe validation/tablas/fase5b_analisis_*.csv y
# validation/figuras/fase5b_cobertura.png. No reajusta nada.
# Uso (desde funresMech_dev):  source("validation/R/16_analisis_fase5b.R")
# ============================================================================

suppressMessages(library(funresMech))
TAB <- "validation/tablas"; RES <- "validation/results/fase5b"; RES5 <- "validation/results/fase5"
FIG <- "validation/figuras"; dir.create(FIG, showWarnings = FALSE)
lee <- function(x) read.csv(file.path(TAB, paste0(x, ".csv")), stringsAsFactors = FALSE)
escribe <- function(x, nombre) write.csv(x, file.path(TAB, paste0("fase5b_analisis_", nombre, ".csv")),
                                         row.names = FALSE)
clave <- function(x) paste(x$id, x$dato)
wilson <- function(x) {
  n <- length(x); p <- mean(x); z <- 1.96
  c(p = p, inf = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    sup = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), n = n)
}
auc <- function(x, y) { ok <- !is.na(x); x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

G <- lee("fase5b_G_simulacion"); E <- lee("fase5b_E_chequeos")
A <- lee("fase5_A_simulacion")
G$lado_fallo <- ifelse(G$cubre, "cubre", ifelse(!is.na(G$z_low) & G$z < G$z_low, "z bajo el IC", "z sobre el IC"))
G$ancho <- G$z_high - G$z_low

# ---- 1. Verificación ---------------------------------------------------------
v <- list()
v$trabajos <- table(G$brazo)
v$rds <- length(list.files(RES, "\\.rds$"))
v$ic_contiene_zhat <- with(G, all((is.na(z_low) | z_low <= z_hat + 1e-9) & (is.na(z_high) | z_hat <= z_high + 1e-9)))
v$fallo_opt <- tapply(G$fallo_opt, G$brazo, sum)
v$puntos_con_piso <- tapply(G$puntos_piso > 0, G$brazo, sum)
v$reintentos <- tapply(G$reintentos_piso, G$brazo, sum)
# G1 usa los mismos datos que A (Fase 5)
mismos <- vapply(list.files(RES, "^G1__.*\\.rds$", full.names = TRUE), function(f) {
  o <- readRDS(f)
  f5 <- file.path(RES5, sprintf("A__%s__d%d__.rds", o$trabajo$id, o$trabajo$dato))
  if (!file.exists(f5)) return(NA)
  identical(o$fit$data, readRDS(f5)$fit$data)
}, logical(1))
v$G1_mismos_datos_que_A <- c(identicos = sum(mismos, na.rm = TRUE), total = length(mismos))
# G2 usa los mismos datos que G1 (d1-d4)
G2d <- vapply(list.files(RES, "^G2__.*\\.rds$", full.names = TRUE), function(f) {
  o <- readRDS(f)
  identical(o$fit$data, readRDS(sub("G2__", "G1__", sub("standard", "exploratory", f)))$fit$data)
}, logical(1))
v$G2_mismos_datos_que_G1 <- c(identicos = sum(G2d), total = length(G2d))
print(v)
escribe(data.frame(chequeo = c("rds", "ic_contiene_zhat", "G1_mismos_datos_que_A", "G2_mismos_datos_que_G1",
                               "fallo_opt_total", "ajustes_con_piso", "reintentos_piso"),
                   valor = c(v$rds, v$ic_contiene_zhat,
                             paste(v$G1_mismos_datos_que_A, collapse = "/"),
                             paste(v$G2_mismos_datos_que_G1, collapse = "/"),
                             sum(v$fallo_opt), sum(v$puntos_con_piso), sum(v$reintentos))), "verificacion")

# ---- 2. Cobertura ------------------------------------------------------------
glob <- do.call(rbind, lapply(split(G, G$brazo), function(q) {
  w <- wilson(q$cubre); data.frame(brazo = q$brazo[1], n = nrow(q), cobertura = w[["p"]],
                                   inf = w[["inf"]], sup = w[["sup"]])
}))
glob <- rbind(glob, { w <- wilson(G$cubre); data.frame(brazo = "todos", n = nrow(G), cobertura = w[["p"]],
                                                      inf = w[["inf"]], sup = w[["sup"]]) })
escribe(glob, "cobertura_global")
porz <- do.call(rbind, lapply(split(G, list(G$brazo, G$z), drop = TRUE), function(q) {
  w <- wilson(q$cubre)
  data.frame(brazo = q$brazo[1], z = q$z[1], n = nrow(q), cobertura = w[["p"]], inf = w[["inf"]],
             sup = w[["sup"]], z_bajo_IC = mean(q$lado_fallo == "z bajo el IC"),
             z_sobre_IC = mean(q$lado_fallo == "z sobre el IC"), ancho_med = median(q$ancho, na.rm = TRUE))
}))
escribe(porz, "cobertura_por_z")
esc <- do.call(rbind, lapply(split(G, list(G$brazo, G$diseno, G$z, G$k), drop = TRUE), function(q)
  data.frame(brazo = q$brazo[1], diseno = q$diseno[1], z = q$z[1], k = q$k[1], n = nrow(q),
             cobertura = mean(q$cubre), z_mediana = median(q$z_hat), ancho_med = median(q$ancho, na.rm = TRUE),
             z_bajo_IC = sum(q$lado_fallo == "z bajo el IC"), z_sobre_IC = sum(q$lado_fallo == "z sobre el IC"))))
escribe(esc, "cobertura_por_escenario")
# factores: k y diseño, todos los brazos con exploratory (G1 + G3)
GE <- G[G$variante == "exploratory", ]
fact <- rbind(
  do.call(rbind, lapply(split(GE, GE$k), function(q) { w <- wilson(q$cubre)
    data.frame(factor = "k", nivel = q$k[1], n = nrow(q), cobertura = w[["p"]], inf = w[["inf"]], sup = w[["sup"]]) })),
  do.call(rbind, lapply(split(GE, GE$diseno), function(q) { w <- wilson(q$cubre)
    data.frame(factor = "diseno", nivel = q$diseno[1], n = nrow(q), cobertura = w[["p"]], inf = w[["inf"]], sup = w[["sup"]]) })))
escribe(fact, "cobertura_por_factor")
print(glob, digits = 3); print(porz, digits = 3); print(fact, digits = 3)

# ---- 3. Comparaciones pareadas -----------------------------------------------
A5 <- A[A$bloque == "nucleo", ]
g1 <- G[G$brazo == "G1", ]
m1 <- merge(g1, A5[, c("id", "dato", "cubre", "z_hat", "z_low", "z_high", "n_puntos", "segundos")],
            by = c("id", "dato"), suffixes = c("", "_f5"))
m1$ancho_f5 <- m1$z_high_f5 - m1$z_low_f5
g2 <- G[G$brazo == "G2", ]
m2 <- merge(g2, g1[, c("id", "dato", "cubre", "z_hat", "ancho", "nll_prec_hat", "segundos")],
            by = c("id", "dato"), suffixes = c("", "_g1"))
mcnemar_p <- function(a, b) { n01 <- sum(!a & b); n10 <- sum(a & !b)
  if (n01 + n10 == 0) return(1); stats::binom.test(min(n01, n10), n01 + n10)$p.value }
par <- data.frame(
  comparacion = c("G1 (1.2.0rc2) vs Fase 5 (mismos datos)", "G2 standard vs G1 exploratory (mismos datos)"),
  n = c(nrow(m1), nrow(m2)),
  cobertura_ref = c(mean(m1$cubre_f5), mean(m2$cubre_g1)), cobertura_nueva = c(mean(m1$cubre), mean(m2$cubre)),
  solo_ref_cubre = c(sum(m1$cubre_f5 & !m1$cubre), sum(m2$cubre_g1 & !m2$cubre)),
  solo_nueva_cubre = c(sum(!m1$cubre_f5 & m1$cubre), sum(!m2$cubre_g1 & m2$cubre)),
  p_mcnemar = c(mcnemar_p(m1$cubre_f5, m1$cubre), mcnemar_p(m2$cubre_g1, m2$cubre)),
  mismo_zhat = c(mean(abs(m1$z_hat - m1$z_hat_f5) < 1e-8), mean(abs(m2$z_hat - m2$z_hat_g1) < 1e-8)),
  ancho_med_ref = c(median(m1$ancho_f5, na.rm = TRUE), median(m2$ancho_g1, na.rm = TRUE)),
  ancho_med_nueva = c(median(m1$ancho, na.rm = TRUE), median(m2$ancho, na.rm = TRUE)),
  horas_ref = c(mean(m1$segundos_f5), mean(m2$segundos_g1)) / 3600,
  horas_nueva = c(mean(m1$segundos), mean(m2$segundos)) / 3600)
par$nll_mejor_g2 <- c(NA, mean(m2$nll_prec_hat < m2$nll_prec_hat_g1 - 0.5))
escribe(par, "pareadas")
print(par, digits = 3)

# ---- 4. Error tipo I, potencia, sesgo ----------------------------------------
G$aic2 <- G$dAIC_z1 > 2
pot <- do.call(rbind, lapply(split(G, list(G$brazo, G$z), drop = TRUE), function(q) {
  w <- wilson(q$excluye_1)
  data.frame(brazo = q$brazo[1], z = q$z[1], n = nrow(q),
             medida = if (q$z[1] == 1) "error tipo I" else "potencia",
             IC_excluye_1 = w[["p"]], inf = w[["inf"]], sup = w[["sup"]], dAIC_mayor_2 = mean(q$aic2),
             sesgo_z_mediana = median(q$z_hat - q$z), h_rel_mediana = median(q$h_hat / (q$hT * 24)),
             k_rel_mediana = median(q$k_hat / q$k), abiertos = mean(q$abierto))
}))
escribe(pot, "tipoI_potencia_sesgo")
print(pot, digits = 3)

# ---- 5. Calibración de los chequeos (H18) -------------------------------------
E$bueno <- E$cubre & !E$fallo_opt & abs(E$err_z) <= 1
th <- funresMech:::.rel_thresholds()
cal <- do.call(rbind, lapply(split(E, list(E$variante, E$check), drop = TRUE), function(q)
  data.frame(preset = q$variante[1], check = q$check[1], n = nrow(q), malos = sum(!q$bueno),
             alerta = mean(q$status != "ok"), falla = mean(q$status == "fail"),
             alerta_buenos = mean(q$status[q$bueno] != "ok"), alerta_malos = mean(q$status[!q$bueno] != "ok"),
             auc_malo = auc(q$value, !q$bueno), auc_no_cubre = auc(q$value, !q$cubre),
             p50_buenos = unname(stats::quantile(q$value[q$bueno], 0.5, na.rm = TRUE)),
             p90_buenos = unname(stats::quantile(q$value[q$bueno], 0.9, na.rm = TRUE)),
             p95_buenos = unname(stats::quantile(q$value[q$bueno], 0.95, na.rm = TRUE)))))
escribe(cal, "calibracion_chequeos")
print(cal[cal$preset == "standard", ], digits = 3)
cat("\nUmbrales actuales:\n"); str(th)

# ---- Figura -------------------------------------------------------------------
png(file.path(FIG, "fase5b_cobertura.png"), 1800, 900, res = 200)
op <- par(mar = c(4.5, 4.5, 2, 1))
cols <- c(G1 = "#1b6ca8", G2 = "#d9822b", G3 = "#3a9a5b")
x0 <- c(1, 1.75, 2.75); off <- c(G1 = -0.08, G2 = 0, G3 = 0.08)
plot(NA, xlim = c(0.7, 3.05), ylim = c(0.5, 1), xaxt = "n", xlab = "True z", ylab = "Coverage of the 95% interval",
     main = "Phase 5b: coverage by arm")
axis(1, x0, x0); abline(h = 0.95, lty = 2)
for (b in names(cols)) { q <- porz[porz$brazo == b, ]
  points(q$z + off[b], q$cobertura, pch = 19, col = cols[b])
  arrows(q$z + off[b], q$inf, q$z + off[b], q$sup, angle = 90, code = 3, length = 0.03, col = cols[b]) }
legend("bottomright", c("G1 core, exploratory", "G2 core, standard", "G3 13 densities, exploratory"),
       col = cols, pch = 19, cex = 0.75, bty = "n")
par(op); dev.off()
cat("\nListo: tablas fase5b_analisis_*.csv y figura", file.path(FIG, "fase5b_cobertura.png"), "\n")
