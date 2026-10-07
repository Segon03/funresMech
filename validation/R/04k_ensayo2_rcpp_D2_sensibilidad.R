# ============================================================
# R/04k_ensayo2_rcpp_D2_sensibilidad.R
#
# Prueba de la hipotesis sobre D2 (Ensayo 2): la replica con 2 de 10
# hospedadores parasitados en N = 10 obliga a un k (forma de la Gamma
# del tiempo de busqueda) muy chico (~0.05, contra ~0.8 en D3). Con k
# chico la probabilidad de encuentro casi no depende de la tasa a*N^z,
# z se compensa con a y k, y el perfil de z queda plano por arriba
# (IC de D2 en escala log: [1.1, >=20], R/04j).
#
# Prediccion: SIN esa replica, k sube hacia valores como los de D3 y el
# perfil de z se cierra (cruza NLL_min + 1.92 en algun z razonable).
#
# Variantes de D2 (escala log, misma configuracion que 04h/04j):
#   sin_2de10   D2 sin la replica (N = 10, parasitados = 2). Perfil
#               completo: 19 valores de z entre 0.75 y 20.
#   control     D2 sin OTRA replica de N = 10 (parasitados = 10), para
#               separar el efecto de "quitar una observacion cualquiera"
#               del de quitar la replica atipica. Perfil reducido: z =
#               2.5 (zona del minimo), 6, 10 y 20.
# El perfil de D2 completo ya existe (04h + 04j) y se usa como referencia.
#
# Costo: 23 tareas (19 + 4) -> una sola tanda con 23 workers, ~2-2.5 h
# de pared (con la configuracion completa de 04h: n_sim 10000, 200
# generaciones, NP 40).
#
# Al terminar (o si ya no quedan tareas) el script resume solo:
#   tablas/ensayo2_rcpp_D2_sensibilidad.csv            una fila por variante
#   tablas/ensayo2_rcpp_D2_sensibilidad_perfiles.csv   todos los puntos
#   figuras/fig_ensayo2_rcpp_D2_sensibilidad.png       perfiles superpuestos
#
# Robustez: igual que 02c (correccion del 04-oct-2026): cache de Rcpp por
# proceso (fun/motor_rcpp.R), una falla en una tarea no corta la corrida,
# escritura atomica de los .rds y recalculo de archivos ilegibles.
#
# REQUISITOS: fun/motor_okuyama.cpp, fun/motor_rcpp.R (version del
# 04-oct-2026), fun/profile_utils.R, data_clean/D2.rds y, para la
# referencia, results/ensayo2_rcpp/log/ y results/ensayo2_rcpp/log_ext/.
# Correr desde la carpeta del proyecto (.Rproj), idealmente como
# Background Job:
#   rstudioapi::jobRunScript("R/04k_ensayo2_rcpp_D2_sensibilidad.R", workingDir = getwd())
# ============================================================

library(future)
library(future.apply)
source(file.path("fun", "motor_rcpp.R"))
source(file.path("fun", "profile_utils.R"))

# ---- Configuracion ----------------------------------------------------------
MODO_PRUEBA <- FALSE   # TRUE: corrida de humo de ~1 min (config minima, pocos z;
                       # resultados en results/ensayo2_rcpp_PRUEBA/sens_D2/)

T_exp    <- 24
N_SIM    <- 10000      # = 04h/04j
ITERMAX  <- 200
NP       <- 40
N_REEVAL <- 5
VARIANTES <- list(
  sin_2de10 = list(quitar = c(dens = 10, par = 2),
                   z = c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 15, 20)),
  control   = list(quitar = c(dens = 10, par = 10),
                   z = c(2.5, 6, 10, 20))
)
N_WORKERS <- max(1, parallel::detectCores() - 1)

if (MODO_PRUEBA) {
  N_SIM <- 300; ITERMAX <- 3; N_REEVAL <- 2
  VARIANTES$sin_2de10$z <- c(2.5, 6)
  VARIANTES$control$z   <- c(2.5)
}
dir_base <- file.path("results", if (MODO_PRUEBA) "ensayo2_rcpp_PRUEBA" else "ensayo2_rcpp", "sens_D2")

# ---- Datos --------------------------------------------------------------------
f_d2 <- file.path("data_clean", "D2.rds")
if (!file.exists(f_d2)) stop("No se encontro ", f_d2)
D2 <- as.data.frame(readRDS(f_d2))

quitar_una <- function(d, dens, par) {
  i <- which(d$dens == dens & d$par == par)
  if (!length(i)) stop(sprintf("No hay ninguna fila con dens = %g y par = %g en D2", dens, par))
  d[-i[1], ]   # quita UNA sola fila
}
datasets <- lapply(VARIANTES, function(v) quitar_una(D2, v$quitar[["dens"]], v$quitar[["par"]]))
for (nm in names(datasets))
  cat(sprintf("Variante %-10s: D2 sin una replica (N = %g, parasitados = %g) -> %d ensayos (D2 completo: %d)\n",
              nm, VARIANTES[[nm]]$quitar[["dens"]], VARIANTES[[nm]]$quitar[["par"]],
              nrow(datasets[[nm]]), nrow(D2)))

# ---- Tareas -------------------------------------------------------------------
tareas <- list()
for (nm in names(VARIANTES)) {
  dir.create(file.path(dir_base, nm), showWarnings = FALSE, recursive = TRUE)
  for (z in VARIANTES[[nm]]$z) {
    f <- file.path(dir_base, nm, sprintf("perfil_D2_z%.2f.rds", z))
    if (file.exists(f)) {
      ok <- tryCatch({ r <- readRDS(f); is.list(r) && is.finite(r$nll) }, error = function(e) FALSE)
      if (isTRUE(ok)) next
      cat("  Archivo ilegible, se recalcula:", f, "\n")
      file.rename(f, paste0(f, ".corrupto"))
    }
    tareas[[length(tareas) + 1]] <- list(variante = nm, z = z, archivo = f)
  }
}
n_total <- sum(sapply(VARIANTES, function(v) length(v$z)))
cat(sprintf("\nTareas: %d de %d pendientes\n", length(tareas), n_total))

if (length(tareas)) {
  cat("Compilando/cargando el motor C++ ...\n")
  cargar_motor_rcpp()
  n_workers <- max(1, min(N_WORKERS, length(tareas)))
  cat(sprintf("Con %d workers: ~%s h de pared (estimacion con los tiempos de 04h)\n\n",
              n_workers, if (MODO_PRUEBA) "0" else
                sprintf("%.0f-%.0f", 2 * ceiling(length(tareas) / n_workers), 2.5 * ceiling(length(tareas) / n_workers))))
  plan(multisession, workers = n_workers)

  correr_tarea <- function(tarea) {
    if (file.exists(tarea$archivo)) return(sprintf("%s z=%.2f: ya existe", tarea$variante, tarea$z))
    un_intento <- function() {
      source(file.path("fun", "motor_rcpp.R"))
      cargar_motor_rcpp()
      res <- ajustar_z_fijo(tarea$z, datasets[[tarea$variante]], T = CFG$T_exp, nsim = CFG$N_SIM,
                            escala = "log", itermax = CFG$ITERMAX, NP = CFG$NP,
                            n_reeval = CFG$N_REEVAL)
      res$dataset <- paste0("D2_", tarea$variante)
      tmp <- paste0(tarea$archivo, ".tmp", Sys.getpid())
      saveRDS(res, tmp)
      if (!file.rename(tmp, tarea$archivo)) stop("no se pudo renombrar ", tmp)
      res
    }
    res <- tryCatch(un_intento(), error = function(e1) {
      Sys.sleep(5)
      tryCatch(un_intento(), error = function(e2) e2)
    })
    if (inherits(res, "error"))
      return(sprintf("ERROR %s z=%.2f: %s", tarea$variante, tarea$z, conditionMessage(res)))
    sprintf("%-9s z=%5.2f: nll=%.3f reeval=%.3f (a=%.3g h=%.3f k=%.3g s=%.3f) %.1f min%s",
            tarea$variante, tarea$z, res$nll, res$nll_reeval, res$par[["a"]], res$par[["h"]],
            res$par[["k"]], res$par[["s"]], res$segundos / 60,
            if (any(res$en_cota)) paste0("  [en cota: ", paste(names(res$en_cota)[res$en_cota], collapse = ","), "]") else "")
  }

  CFG <- list(T_exp = T_exp, N_SIM = N_SIM, ITERMAX = ITERMAX, NP = NP, N_REEVAL = N_REEVAL)
  t0 <- Sys.time()
  mensajes <- future_lapply(tareas, correr_tarea, future.seed = TRUE, future.scheduling = Inf,
                            future.globals = list(datasets = datasets, CFG = CFG),
                            future.packages = c("Rcpp", "DEoptim"))
  invisible(lapply(mensajes, function(m) cat(m, "\n")))
  n_err <- sum(startsWith(unlist(mensajes), "ERROR"))
  if (n_err > 0) cat(sprintf("\nATENCION: %d tareas fallaron. Volver a correr este script las retoma.\n", n_err))
  cat(sprintf("\nAjustes listos en %.1f horas.\n", as.numeric(difftime(Sys.time(), t0, units = "hours"))))
  plan(sequential)
}

# ---- Resumen ------------------------------------------------------------------
leer_dir <- function(d, variante) {
  fs <- list.files(d, pattern = "^perfil_D2_.*\\.rds$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, function(f) {
    r <- readRDS(f)
    data.frame(variante = variante, z = round(r$z, 2), nll = r$nll, nll_reeval = r$nll_reeval,
               a = r$par[["a"]], h = r$par[["h"]], k = r$par[["k"]], s = r$par[["s"]],
               en_cota = paste(names(r$en_cota)[r$en_cota], collapse = ","),
               n_obs = NA_integer_, stringsAsFactors = FALSE)
  }))
}
dir_ref <- file.path("results", "ensayo2_rcpp")
perf <- rbind(
  leer_dir(file.path(dir_ref, "log"), "completo"),
  leer_dir(file.path(dir_ref, "log_ext"), "completo"),
  leer_dir(file.path(dir_base, "sin_2de10"), "sin_2de10"),
  leer_dir(file.path(dir_base, "control"), "control"))
if (is.null(perf) || !nrow(perf)) stop("No hay perfiles para resumir.")
perf$n_obs <- ifelse(perf$variante == "completo", nrow(D2), nrow(D2) - 1L)
perf <- perf[order(perf$variante, perf$z), ]
perf <- perf[!duplicated(perf[, c("variante", "z")]), ]

res <- list()
for (v in c("completo", "sin_2de10", "control")) {
  p <- perf[perf$variante == v, ]
  if (nrow(p) < 1) next
  i <- which.min(p$nll); d <- p$nll - p$nll[i]
  ci <- if (nrow(p) >= 3) ci_from_profile_okuyama(p$z, p$nll) else NULL
  en <- function(z0) { j <- which(abs(p$z - z0) < 1e-9); if (length(j)) d[j] else NA_real_ }
  res[[v]] <- data.frame(
    variante = v, n_obs = p$n_obs[1], n_puntos = nrow(p), z_max = max(p$z),
    z_hat = p$z[i], nll_min = p$nll[i], aic = 2 * p$nll[i] + 10,
    ic_low  = if (is.null(ci)) NA else if (ci$low_censurado) NA else ci$z_low,
    ic_high = if (is.null(ci)) NA else if (ci$high_censurado) NA else ci$z_high,
    ic_low_censurado  = if (is.null(ci)) NA else ci$low_censurado,
    ic_high_censurado = if (is.null(ci)) NA else ci$high_censurado,
    a_zhat = p$a[i], h_zhat = p$h[i], k_zhat = p$k[i], s_zhat = p$s[i],
    dnll_z6 = en(6), dnll_z10 = en(10), dnll_z20 = en(20),
    k_z6 = if (any(p$z == 6)) p$k[p$z == 6] else NA, k_z20 = if (any(p$z == 20)) p$k[p$z == 20] else NA,
    puntos_en_cota = sum(nzchar(p$en_cota)), stringsAsFactors = FALSE)
}
tabla <- do.call(rbind, res)
dir.create("tablas", showWarnings = FALSE); dir.create("figuras", showWarnings = FALSE)
pref <- if (MODO_PRUEBA) "PRUEBA_" else ""
write.csv(tabla, file.path("tablas", paste0(pref, "ensayo2_rcpp_D2_sensibilidad.csv")), row.names = FALSE)
write.csv(perf,  file.path("tablas", paste0(pref, "ensayo2_rcpp_D2_sensibilidad_perfiles.csv")), row.names = FALSE)

cat("\n=== D2: sensibilidad a la replica 2/10 en N = 10 (escala log) ===\n")
cat("(dNLL = NLL - NLL_min de la variante; el IC se cierra si dNLL supera 1.92)\n\n")
num <- sapply(tabla, is.numeric); t2 <- tabla; t2[num] <- lapply(t2[num], function(x) signif(x, 3))
print(t2[, c("variante", "n_obs", "n_puntos", "z_hat", "ic_low", "ic_high", "ic_high_censurado",
             "k_zhat", "k_z6", "k_z20", "dnll_z6", "dnll_z10", "dnll_z20", "puntos_en_cota")], row.names = FALSE)
cat("\nNota: 'control' tiene solo 4 puntos (z = 2.5, 6, 10, 20); su z_hat y su IC no son\n",
    "comparables, sirve para ver si quitar una replica CUALQUIERA de N = 10 cambia la\n",
    "forma del perfil por arriba (dNLL en z = 6, 10, 20).\n")

png(file.path("figuras", paste0(pref, "fig_ensayo2_rcpp_D2_sensibilidad.png")), width = 1500, height = 900, res = 170)
par(mar = c(4.2, 4.5, 2.5, 1))
cols <- c(completo = "#1f5fa8", sin_2de10 = "#d9822b", control = "grey45")
ltys <- c(completo = 1, sin_2de10 = 1, control = 2)
ymax <- max(4, min(15, max(unlist(lapply(split(perf, perf$variante), function(p) p$nll - min(p$nll))))))
plot(NA, xlim = c(0.5, 20), ylim = c(0, ymax), xlab = "z", ylab = expression(NLL - NLL[min]),
     main = "D2, escala log: con y sin la replica 2/10 en N = 10", cex.main = 0.95)
abline(h = UMBRAL_IC95, col = "firebrick", lty = 2)
for (v in names(cols)) {
  p <- perf[perf$variante == v, ]
  if (nrow(p)) lines(p$z, p$nll - min(p$nll), col = cols[v], lty = ltys[v], lwd = 2, type = "o", pch = 16, cex = 0.6)
}
legend("topright", c("D2 completo (04h + 04j)", "D2 sin la replica 2/10", "control: D2 sin una replica 10/10",
                     "umbral IC95 (+1.92)"), col = c(cols, "firebrick"), lty = c(ltys, 2), lwd = c(2, 2, 2, 1),
       bty = "o", bg = "white", box.col = "grey80", cex = 0.75)
dev.off()
cat(sprintf("\nArchivos: tablas/%sensayo2_rcpp_D2_sensibilidad.csv, tablas/%sensayo2_rcpp_D2_sensibilidad_perfiles.csv,\n  figuras/%sfig_ensayo2_rcpp_D2_sensibilidad.png\n",
            pref, pref, pref))
