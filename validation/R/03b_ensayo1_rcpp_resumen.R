# ============================================================
# R/03b_ensayo1_rcpp_resumen.R
#
# Resume el Ensayo 1 (definitivo, motor Rcpp + escala log) a partir de
# los perfiles de R/02c_ensayo1_rcpp.R. Rapido (segundos); se puede
# correr con la corrida a medio terminar (usa los puntos que existan y
# avisa cuantos faltan).
#
# Por dataset (escenario x replica), para la NLL de DEoptim y para la
# NLL re-evaluada:
#   - z_hat = minimo del perfil; a, h, k, s en z_hat
#   - IC95% de z (umbral +1.92) interpolado y sobre la grilla, con
#     marca de censura (el perfil no cruza el umbral dentro de la grilla)
#   - cobertura: z verdadero dentro del IC (un lado censurado se toma
#     como abierto: <= 0.25 o >= 5.0)
# Por escenario:
#   - media, sesgo, RMSE y mediana de z_hat; error relativo mediano de
#     a, h, k, s (estimado/verdadero)
#   - cobertura del IC95% de z (con IC binomial exacto, por el n chico)
#   - % de IC con limite superior / inferior no acotado; ancho mediano
#
# Salidas:
#   tablas/ensayo1_rcpp_ajustes.csv    una fila por dataset x tipo de NLL
#   tablas/ensayo1_rcpp_resumen.csv    una fila por escenario x tipo de NLL
#   tablas/ensayo1_rcpp_perfiles.csv   todos los puntos de los perfiles
#   figuras/fig_ensayo1_rcpp_perfiles.png    perfiles por escenario
#   figuras/fig_ensayo1_rcpp_parametros.png  estimado/verdadero por parametro
# ============================================================

source(file.path("fun", "profile_utils.R"))

MODO_PRUEBA <- FALSE   # TRUE para leer la corrida de humo de 02c
dir_res <- file.path("results", if (MODO_PRUEBA) "ensayo1_rcpp_PRUEBA" else "ensayo1_rcpp")
if (!dir.exists(file.path(dir_res, "datos"))) stop("No hay resultados en ", dir_res, ". Correr 02c primero.")

# ---- Verdad y perfiles --------------------------------------------------------
f_datos <- list.files(file.path(dir_res, "datos"), pattern = "\\.rds$", full.names = TRUE)
info <- do.call(rbind, lapply(f_datos, function(f) {
  d <- readRDS(f)
  data.frame(id = sub("\\.rds$", "", basename(f)), escenario = d$escenario, replica = d$replica,
             a_true = d$verdad[["a"]], h_true = d$verdad[["h"]], z_true = d$verdad[["z"]],
             k_true = d$verdad[["k"]], s_true = d$verdad[["s"]], stringsAsFactors = FALSE)
}))

filas <- list()
for (id in info$id) for (f in list.files(file.path(dir_res, id), pattern = "^perfil_.*\\.rds$", full.names = TRUE)) {
  r <- readRDS(f)
  filas[[length(filas) + 1]] <- data.frame(
    id = id, z = round(r$z, 2), nll = r$nll, nll_reeval = r$nll_reeval,
    nll_reeval_se = r$nll_reeval_se,
    a = r$par[["a"]], h = r$par[["h"]], k = r$par[["k"]], s = r$par[["s"]],
    en_cota = paste(names(r$en_cota)[r$en_cota], collapse = ","),
    minutos = r$segundos / 60, nsim = r$config$nsim, itermax = r$config$itermax,
    stringsAsFactors = FALSE)
}
if (!length(filas)) stop("Todavia no hay puntos de perfil en ", dir_res)
perf <- merge(do.call(rbind, filas), info[, c("id", "escenario", "replica")], by = "id")
perf <- perf[order(perf$escenario, perf$replica, perf$z), ]
dir.create("tablas", showWarnings = FALSE); dir.create("figuras", showWarnings = FALSE)
write.csv(perf, file.path("tablas", "ensayo1_rcpp_perfiles.csv"), row.names = FALSE)

z_grid_max <- max(table(perf$id))
cat(sprintf("Perfiles leidos: %d datasets, %d puntos en total (%d por perfil completo)\n",
            length(unique(perf$id)), nrow(perf), z_grid_max))
incompletos <- names(which(table(perf$id) < z_grid_max))
if (length(incompletos)) cat(sprintf("  AVISO: %d perfiles incompletos (corrida sin terminar)\n", length(incompletos)))
if (any(nzchar(perf$en_cota)))
  cat(sprintf("  Aviso: %d puntos con algun parametro en la cota de busqueda (%s)\n",
              sum(nzchar(perf$en_cota)), paste(unique(perf$en_cota[nzchar(perf$en_cota)]), collapse = "; ")))

# ---- Un resultado por dataset x tipo de NLL -------------------------------------
res <- list()
for (id in unique(perf$id)) {
  p <- perf[perf$id == id, ]; v0 <- info[info$id == id, ]
  if (nrow(p) < 3) next
  for (tipo in c("nll", "nll_reeval")) {
    v <- p[[tipo]]; i <- which.min(v)
    ci  <- ci_from_profile_okuyama(p$z, v)
    cig <- ci_grilla_okuyama(p$z, v)
    lo <- if (ci$low_censurado) -Inf else ci$z_low
    hi <- if (ci$high_censurado) Inf else ci$z_high
    res[[length(res) + 1]] <- data.frame(
      id = id, escenario = v0$escenario, replica = v0$replica, tipo_nll = tipo,
      n_puntos = nrow(p), z_true = v0$z_true, z_hat = p$z[i], nll_min = v[i],
      ic_low = ci$z_low, ic_high = ci$z_high,
      ic_low_censurado = ci$low_censurado, ic_high_censurado = ci$high_censurado,
      ic_grilla_low = cig$interior_low, ic_grilla_high = cig$interior_high,
      cubre = v0$z_true >= lo && v0$z_true <= hi,
      a = p$a[i], h = p$h[i], k = p$k[i], s = p$s[i],
      a_true = v0$a_true, h_true = v0$h_true, k_true = v0$k_true, s_true = v0$s_true,
      stringsAsFactors = FALSE)
  }
}
aj <- do.call(rbind, res)
write.csv(aj, file.path("tablas", "ensayo1_rcpp_ajustes.csv"), row.names = FALSE)

# ---- Resumen por escenario --------------------------------------------------------
ic_binom <- function(x, n) if (n > 0) stats::binom.test(x, n)$conf.int else c(NA, NA)
resu <- list()
for (esc in sort(unique(aj$escenario))) for (tipo in c("nll", "nll_reeval")) {
  d <- aj[aj$escenario == esc & aj$tipo_nll == tipo, ]
  if (!nrow(d)) next
  zt <- d$z_true[1]; n <- nrow(d); cb <- sum(d$cubre); icb <- ic_binom(cb, n)
  ancho <- ifelse(d$ic_low_censurado | d$ic_high_censurado, NA, d$ic_high - d$ic_low)
  rel <- function(est, tru) median(est / tru)
  resu[[length(resu) + 1]] <- data.frame(
    escenario = esc, tipo_nll = tipo, n = n, z_true = zt,
    z_hat_media = mean(d$z_hat), z_hat_mediana = median(d$z_hat),
    sesgo_z = mean(d$z_hat) - zt, rmse_z = sqrt(mean((d$z_hat - zt)^2)),
    cobertura = cb / n, cobertura_ic95_inf = icb[1], cobertura_ic95_sup = icb[2],
    pct_ic_sup_no_acotado = 100 * mean(d$ic_high_censurado),
    pct_ic_inf_no_acotado = 100 * mean(d$ic_low_censurado),
    ancho_ic_mediano = if (all(is.na(ancho))) NA else median(ancho, na.rm = TRUE),
    a_rel_mediana = rel(d$a, d$a_true), h_rel_mediana = rel(d$h, d$h_true),
    k_rel_mediana = rel(d$k, d$k_true), s_rel_mediana = rel(d$s, d$s_true),
    stringsAsFactors = FALSE)
}
tabla <- do.call(rbind, resu)
write.csv(tabla, file.path("tablas", "ensayo1_rcpp_resumen.csv"), row.names = FALSE)

cat("\n=== Ensayo 1: recuperacion de parametros (motor Rcpp, escala log) ===\n")
cat("(rel = mediana de estimado/verdadero; 1 = sin sesgo)\n\n")
print(transform(tabla[, c("escenario", "tipo_nll", "n", "z_true", "z_hat_mediana", "sesgo_z", "rmse_z",
                          "cobertura", "pct_ic_sup_no_acotado", "pct_ic_inf_no_acotado", "ancho_ic_mediano",
                          "a_rel_mediana", "h_rel_mediana", "k_rel_mediana", "s_rel_mediana")],
                sesgo_z = round(sesgo_z, 2), rmse_z = round(rmse_z, 2), cobertura = round(cobertura, 2),
                ancho_ic_mediano = round(ancho_ic_mediano, 2), a_rel_mediana = round(a_rel_mediana, 2),
                h_rel_mediana = round(h_rel_mediana, 2), k_rel_mediana = round(k_rel_mediana, 2),
                s_rel_mediana = round(s_rel_mediana, 2)), row.names = FALSE)

# ---- Figuras ------------------------------------------------------------------------
escs <- sort(unique(perf$escenario))
png(file.path("figuras", "fig_ensayo1_rcpp_perfiles.png"), width = 700 * length(escs), height = 800, res = 170)
op <- par(mfrow = c(1, length(escs)), mar = c(4.2, 4.5, 2.5, 1))
for (esc in escs) {
  ids <- unique(perf$id[perf$escenario == esc]); zt <- info$z_true[info$escenario == esc][1]
  plot(NA, xlim = range(perf$z), ylim = c(0, 8), xlab = "z", ylab = expression(NLL - NLL[min]),
       main = sprintf("Escenario %s (z verdadero = %.2f)", esc, zt), cex.main = 0.95)
  for (id in ids) { p <- perf[perf$id == id, ]
    lines(p$z, p$nll - min(p$nll), col = adjustcolor("#1f5fa8", 0.45), lwd = 1) }
  abline(v = zt, col = "grey20", lty = 2, lwd = 1.5)
  abline(h = UMBRAL_IC95, col = "firebrick", lty = 2)
  legend("top", c("perfil de cada replica", "z verdadero", "umbral IC95 (+1.92)"),
         col = c("#1f5fa8", "grey20", "firebrick"), lty = c(1, 2, 2), bty = "n", cex = 0.75)
}
par(op); dev.off()

d1 <- aj[aj$tipo_nll == "nll", ]
png(file.path("figuras", "fig_ensayo1_rcpp_parametros.png"), width = 1900, height = 800, res = 170)
op <- par(mfrow = c(1, 5), mar = c(4, 4.2, 2.5, 0.8))
for (pp in c("z", "a", "h", "k", "s")) {
  est <- if (pp == "z") d1$z_hat else d1[[pp]]
  r <- est / d1[[paste0(pp, "_true")]]
  boxplot(r ~ d1$escenario, log = if (pp %in% c("a", "k")) "y" else "", xlab = "escenario",
          ylab = "estimado / verdadero", main = pp, col = "#dbe6f4", border = "#1f5fa8", outline = TRUE)
  abline(h = 1, col = "firebrick", lty = 2)
}
par(op); dev.off()

cat("\nArchivos: tablas/ensayo1_rcpp_ajustes.csv, tablas/ensayo1_rcpp_resumen.csv,",
    "tablas/ensayo1_rcpp_perfiles.csv, figuras/fig_ensayo1_rcpp_perfiles.png,",
    "figuras/fig_ensayo1_rcpp_parametros.png\n")
