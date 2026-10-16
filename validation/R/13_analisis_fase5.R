# ============================================================================
# 13_analisis_fase5.R
#
# Verificación e interpretación de la Fase 5 (12_validacion_ampliada.R).
# Lee validation/tablas/fase5_*.csv y validation/results/fase5/*.rds y escribe
#   validation/tablas/fase5_analisis_*.csv
#   validation/figuras/fase5_*.png
# No vuelve a ajustar nada: solo lee los resultados guardados.
# Uso (desde funresMech_dev):  source("validation/R/13_analisis_fase5.R")
# ============================================================================

suppressMessages(library(funresMech))
TAB <- "validation/tablas"; RES <- "validation/results/fase5"; FIG <- "validation/figuras"
dir.create(FIG, showWarnings = FALSE)
lee <- function(x) read.csv(file.path(TAB, paste0("fase5_", x, ".csv")), stringsAsFactors = FALSE)
escribe <- function(x, nombre) write.csv(x, file.path(TAB, paste0("fase5_analisis_", nombre, ".csv")),
                                         row.names = FALSE)
clave <- function(x) paste(x$id, x$dato)
prop_ic <- function(x) {                     # proporción con IC de Wilson al 95%
  n <- length(x); p <- mean(x); z <- 1.96
  c(p = p, inf = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    sup = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), n = n)
}

A  <- lee("A_simulacion"); B <- lee("B_misespecificacion"); C <- lee("C_presupuesto")
D  <- lee("D_110_vs_120_simulados"); DR <- lee("D_real"); E <- lee("E_chequeos")
FF <- lee("F_D2_D3_final"); FP <- lee("F_perfiles")

# ---- Paso 1: verificación ----------------------------------------------------
v <- list()
v$n_rds <- length(list.files(RES, "\\.rds$")) - length(list.files(RES, "_ajuste\\.rds$"))
v$conteos <- c(A = nrow(A), B = nrow(B), C = nrow(C), D = nrow(D), Dreal = nrow(DR),
               E = nrow(E), F = nrow(FF))
Ce <- C[C$variante == "exploratory", ]
v$C_igual_A <- isTRUE(all.equal(Ce$nll_prec_hat, A$nll_prec_hat[match(clave(Ce), clave(A))]))
D12 <- D[D$version == "1.2.0", ]
v$D_igual_A <- isTRUE(all.equal(D12$z_hat, A$z_hat[match(clave(D12), clave(A))]))
v$ic_contiene_zhat <- with(A, all((is.na(z_low) | z_low <= z_hat) & (is.na(z_high) | z_hat <= z_high)))
v$deficit_mayor_1 <- mean(A$deficit > 1)
v$F_sin_atipico <- FF$ensayos[FF$datos == "D2_sin_atipico"] == 72 &&
  FF$filas_quitadas[FF$datos == "D2_sin_atipico"] == "13"
v$F_min_perfil_es_zhat <- all(vapply(FF$datos, function(d) {
  q <- FP[FP$datos == d, ]; q$z[which.min(q$nll)] == FF$z_hat[FF$datos == d] }, logical(1)))
print(v)

# Perfiles de A (desde los .rds): saltos por piso de probabilidad
fa <- list.files(RES, "^A__.*\\.rds$", full.names = TRUE)
PA <- do.call(rbind, lapply(fa, function(f) {
  o <- readRDS(f); p <- o$fit$profile
  data.frame(id = o$trabajo$id, dato = o$trabajo$dato, p[, c("z", "nll", "nll_se", "nll_opt", "n_floor")])
}))
PA$salto <- PA$nll - PA$nll_opt
PA$z_true <- as.numeric(sub(".*_z([0-9.]+)_k.*", "\\1", PA$id))
saltos <- data.frame(puntos = nrow(PA), puntos_con_piso = sum(PA$n_floor > 0),
                     ajustes_con_piso = length(unique(clave(PA)[PA$n_floor > 0])),
                     ajustes_con_piso_en_z1 = length(unique(clave(PA)[PA$n_floor > 0 & PA$z == 1])))
escribe(saltos, "pisos_perfil")

# La NLL precisa en la verdad también puede tocar el piso (déficits muy negativos)
v$verdad_con_piso <- sum(A$nll_prec_true - A$nll > 50)

# ---- Paso 2A: simulación amplia ----------------------------------------------
A$err_z <- A$z_hat - A$z
A$aic2 <- A$dAIC_z1 > 2
A$abierto <- A$abierto_bajo | A$abierto_alto
A$h_rel <- A$h_hat / (A$hT * 24)
resA <- do.call(rbind, lapply(split(A, list(A$bloque, A$diseno, A$z), drop = TRUE), function(q)
  data.frame(bloque = q$bloque[1], diseno = q$diseno[1], z = q$z[1], n = nrow(q),
             cobertura = mean(q$cubre), abiertos = mean(q$abierto),
             z_mediana = median(q$z_hat), err_abs_med = median(abs(q$err_z)),
             excluye_1 = mean(q$excluye_1), aic_mayor_2 = mean(q$aic2),
             h_rel_med = median(q$h_rel), fallo_opt = mean(q$deficit > 1))))
escribe(resA, "A_por_escenario")
globalA <- rbind(
  cbind(medida = "cobertura (todas)", t(prop_ic(A$cubre))),
  cbind(medida = "cobertura z = 1", t(prop_ic(A$cubre[A$z == 1]))),
  cbind(medida = "cobertura z = 1.75", t(prop_ic(A$cubre[A$z == 1.75]))),
  cbind(medida = "cobertura z = 2.75", t(prop_ic(A$cubre[A$z == 2.75]))),
  cbind(medida = "error tipo I (IC excluye 1, z = 1)", t(prop_ic(A$excluye_1[A$z == 1]))),
  cbind(medida = "potencia z = 1.75", t(prop_ic(A$excluye_1[A$z == 1.75]))),
  cbind(medida = "potencia z = 2.75", t(prop_ic(A$excluye_1[A$z == 2.75]))),
  cbind(medida = "IC abierto", t(prop_ic(A$abierto))),
  cbind(medida = "falla de optimización (déficit > 1)", t(prop_ic(A$deficit > 1))))
escribe(as.data.frame(globalA), "A_global")

# Causas de la subcobertura: IC recalculado sobre los mismos perfiles
#   (1) sin los puntos con piso de probabilidad; (2) además con interpolación
#   por spline monótona en vez de lineal.
ci_spline <- function(z, nll) {
  f <- stats::splinefun(z, nll, method = "monoH.FC")
  zz <- seq(min(z), max(z), by = 0.01)
  funresMech:::.profile_ci(zz, f(zz))
}
cubre <- function(ci, zt) (is.na(ci$z_low) || ci$z_low <= zt) && (is.na(ci$z_high) || zt <= ci$z_high)
alt <- do.call(rbind, lapply(split(PA, clave(PA)), function(q) {
  q <- q[order(q$z), ]; zt <- q$z_true[1]; q2 <- q[q$n_floor == 0, ]
  c0 <- funresMech:::.profile_ci(q$z, q$nll); c1 <- funresMech:::.profile_ci(q2$z, q2$nll)
  data.frame(z = zt, actual = cubre(c0, zt), sin_piso = cubre(c1, zt),
             sin_piso_spline = cubre(ci_spline(q2$z, q2$nll), zt))
}))
causas <- aggregate(cbind(actual, sin_piso, sin_piso_spline) ~ z, alt, mean)
causas <- rbind(causas, data.frame(z = NA, t(colMeans(alt[, -1]))))
escribe(causas, "A_cobertura_alternativas")

# ---- Paso 2B: desvíos ---------------------------------------------------------
B$err_z <- B$z_hat - B$z; B$h_rel <- B$h_hat / (B$hT * 24)
resB <- aggregate(cbind(cubre, err_z, excluye_1, h_rel, dAIC_BB) ~ gen + variante + z, B, mean)
escribe(resB, "B_resumen")
w <- merge(B[B$variante == "ajuste_por_defecto", c("id", "dato", "nll")],
           B[B$variante == "ajuste_correcto", c("id", "dato", "nll", "gen")], by = c("id", "dato"))
w$dAIC <- 2 * (w$nll.x - w$nll.y)
selB <- aggregate(dAIC ~ gen, w, function(x) c(mediana = median(x), correcto_mejor = mean(x > 0),
                                               dAIC_mayor_2 = mean(x > 2)))
escribe(do.call(data.frame, selB), "B_seleccion_familia")

# ---- Paso 2C: presupuesto ----------------------------------------------------
C$abierto <- C$abierto_bajo | C$abierto_alto
resC <- aggregate(cbind(exito_05, exito_1, brecha, cubre, abierto, segundos) ~ variante, C, mean)
escribe(resC, "C_resumen")

# ---- Paso 2D: 1.1.0 frente a 1.2.0 -------------------------------------------
D$abierto <- D$abierto_bajo | D$abierto_alto
resD <- aggregate(cbind(cubre, abierto, abs_err_z, excluye_1, brecha, segundos) ~ version, D, mean)
escribe(resD, "D_resumen")
resDR <- aggregate(cbind(z_hat, z_low, nll_prec_hat, segundos) ~ id + variante, DR,
                   function(x) c(media = mean(x), rango = diff(range(x))))
escribe(do.call(data.frame, resDR), "D_real_resumen")

# ---- Paso 2E: chequeos -------------------------------------------------------
A$malo <- !A$cubre | A$deficit > 1 | abs(A$err_z) > 1
auc <- function(x, y) { ok <- !is.na(x); x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
k <- clave(A)
resE <- do.call(rbind, lapply(split(E, E$check), function(q) {
  m <- match(clave(q), k)
  data.frame(check = q$check[1], alertas = mean(q$status != "ok"), fallas = mean(q$status == "fail"),
             auc_no_cubre = auc(q$value, !A$cubre[m]), auc_fallo_opt = auc(q$value, A$deficit[m] > 1),
             auc_malo = auc(q$value, A$malo[m]),
             p90_buenos = unname(stats::quantile(q$value[!A$malo[m]], 0.9, na.rm = TRUE)))
}))
escribe(resE, "E_discriminacion")

# ---- Paso 2F: D2 y D3 --------------------------------------------------------
escribe(FF[, c("datos", "ensayos", "filas_quitadas", "z_hat", "z_low", "z_high", "a_hat", "h_hat",
               "k_hat", "s_hat", "nll", "nll_se", "dAIC_z1", "dAIC_BB", "n_puntos", "fiabilidad")],
        "F_resumen")

# ---- Figuras ------------------------------------------------------------------
png(file.path(FIG, "fase5_cobertura.png"), 1800, 1100, res = 200)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))
cz <- t(sapply(split(A$cubre, A$z), prop_ic))
x <- seq_len(nrow(cz))
plot(x, cz[, "p"], ylim = c(0.5, 1), xaxt = "n", pch = 19, xlab = "True z", ylab = "Coverage of the 95% interval",
     main = "A. Coverage by true z")
axis(1, x, rownames(cz)); arrows(x, cz[, "inf"], x, cz[, "sup"], angle = 90, code = 3, length = 0.04)
abline(h = 0.95, lty = 2)
m <- as.matrix(causas[!is.na(causas$z), -1]); rownames(m) <- causas$z[!is.na(causas$z)]
barplot(t(m), beside = TRUE, ylim = c(0, 1), legend.text = c("current", "without floor points",
        "without floor + spline"), args.legend = list(x = "bottomright", cex = 0.7),
        xlab = "True z", ylab = "Coverage", main = "Recomputed on the same profiles")
abline(h = 0.95, lty = 2); par(op); dev.off()

png(file.path(FIG, "fase5_perfiles_D2_D3.png"), 1800, 700, res = 200)
op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1))
for (d in c("D2", "D2_sin_atipico", "D3")) {
  q <- FP[FP$datos == d & FP$z <= 10, ]
  plot(q$z, q$nll - min(q$nll), type = "b", pch = 19, cex = 0.6, xlab = "z", ylab = "NLL - min",
       main = d, ylim = c(0, 6)); abline(h = 1.92, lty = 2)
}
par(op); dev.off()
cat("\nListo: tablas fase5_analisis_*.csv y figuras en", FIG, "\n")
