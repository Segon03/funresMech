# ============================================================
# R/06_diagnostico_atipicos.R
#
# Diagnostico de ensayos atipicos/influyentes (fun/diagnostico_datos.R)
# aplicado a:
#   (1) D2 y D3 (Ensayo 2): cribado leave-one-out ANTES del ajuste, con
#       la Beta-Binomial de media Okuyama (rapido, segundos).
#   (2) Los 60 datasets simulados del Ensayo 1 (sin atipicos por
#       construccion): tasa de falsos positivos del cribado.
#   (3) Residuos cuantilicos con la predictiva SIMULADA (motor C++) en
#       los parametros ajustados de D2 (con y sin la replica 2/10) y D3:
#       muestra el efecto de enmascaramiento (el ajuste con k chico
#       "absorbe" al ensayo atipico).
#
# Resultado de referencia (05-oct-2026):
#   - D2: un solo ensayo marcado, N = 10 con 2 parasitados (p LOO =
#     5e-7, umbral Bonferroni 6.8e-4); sin el, la sobredispersion cae
#     97%. D3: ninguno.
#   - Ensayo 1: 0/20 datasets con falsos positivos en A y en B; en C
#     (k chico, colas pesadas por construccion) 4/20 datasets con un
#     ensayo marcado (p entre 2e-6 y 1e-4; D2 es mas extremo: 5e-7).
#   - Residuos: con el ajuste de D2 completo el ensayo 2/10 tiene
#     P(Y<=2) = 0.014 (residuo -2.3, no se ve raro: enmascaramiento);
#     con el ajuste sin el, P = 5e-6 (residuo -4.5).
#
# Duracion: ~10-15 min con varios nucleos (el cribado del Ensayo 1 hace
# 60 x 65 ajustes Beta-Binomiales). Necesita: fun/betabinom_model.R,
# fun/diagnostico_datos.R, fun/motor_rcpp.R (+ Rcpp/Rtools), extraDistr,
# future.apply; data_clean/D2.rds, D3.rds; results/ensayo1_rcpp/datos/.
#
# Salidas:
#   tablas/diagnostico_cribado_D2_D3.csv
#   tablas/diagnostico_cribado_ensayo1.csv
#   tablas/diagnostico_residuos_pit.csv
# ============================================================

library(future)
library(future.apply)
source(file.path("fun", "betabinom_model.R"))
source(file.path("fun", "diagnostico_datos.R"))
source(file.path("fun", "motor_rcpp.R"))
T_exp <- 24
dir.create("tablas", showWarnings = FALSE)

# ---- (1) Cribado de D2 y D3 ------------------------------------------------
cr <- list()
for (lab in c("D2", "D3")) {
  d <- as.data.frame(readRDS(file.path("data_clean", paste0(lab, ".rds"))))
  r <- cribado_atipicos(d, T = T_exp)
  r$dataset <- lab; r$rho_todos <- attr(r, "rho_todos"); r$umbral_p <- attr(r, "umbral_p")
  cr[[lab]] <- r
  cat(sprintf("\n== Cribado %s (n = %d; umbral p = %.1e; rho con todos = %.4f)\n",
              lab, nrow(d), attr(r, "umbral_p"), attr(r, "rho_todos")))
  print(head(r[, c("fila", "dens", "par", "esperado_loo", "p_loo", "cambio_rho", "atipico",
                   "influye_dispersion")], 5), row.names = FALSE, digits = 3)
}
write.csv(do.call(rbind, cr), file.path("tablas", "diagnostico_cribado_D2_D3.csv"), row.names = FALSE)

# ---- (2) Falsos positivos en el Ensayo 1 -------------------------------------
dir_e1 <- file.path("results", "ensayo1_rcpp", "datos")
if (dir.exists(dir_e1)) {
  fs <- list.files(dir_e1, pattern = "\\.rds$", full.names = TRUE)
  plan(multisession, workers = max(1, min(length(fs), parallel::detectCores() - 1)))
  e1 <- future_lapply(fs, function(f) {
    source(file.path("fun", "betabinom_model.R")); source(file.path("fun", "diagnostico_datos.R"))
    x <- readRDS(f); r <- cribado_atipicos(x$datos, T = x$T)
    data.frame(id = sub("\\.rds$", "", basename(f)), escenario = x$escenario,
               n_marcados = sum(r$atipico), n_influye_dispersion = sum(r$influye_dispersion),
               p_min = min(r$p_loo), ensayo_p_min = paste0(r$dens[1], "/", r$par[1]),
               rho = attr(r, "rho_todos"))
  }, future.seed = TRUE)
  plan(sequential)
  e1 <- do.call(rbind, e1)
  write.csv(e1, file.path("tablas", "diagnostico_cribado_ensayo1.csv"), row.names = FALSE)
  cat("\n== Falsos positivos del cribado en el Ensayo 1 (datasets con >= 1 ensayo marcado)\n")
  print(aggregate(cbind(con_marcados = n_marcados > 0) ~ escenario, e1, sum))
  cat(sprintf("p minimo en C: %.1e .. %.1e (D2: %.1e)\n", min(e1$p_min[e1$escenario == "C"]),
              max(e1$p_min[e1$escenario == "C" & e1$n_marcados > 0]), min(cr$D2$p_loo)))
}

# ---- (3) Residuos cuantilicos con la predictiva simulada ----------------------
cargar_motor_rcpp()
leer_par <- function(f) { r <- readRDS(f); list(z = r$z, p = unname(r$par[c("a", "h", "k", "s")])) }
ajustes <- list(
  D2_completo  = c(datos = "D2", archivo = file.path("results", "ensayo2_rcpp", "log", "perfil_D2_z2.80.rds")),
  D2_sin_2de10 = c(datos = "D2", archivo = file.path("results", "ensayo2_rcpp", "sens_D2", "sin_2de10", "perfil_D2_z2.00.rds")),
  D3           = c(datos = "D3", archivo = file.path("results", "ensayo2_rcpp", "log", "perfil_D3_z2.50.rds")))
res_pit <- list()
for (nm in names(ajustes)) {
  a <- ajustes[[nm]]
  if (!file.exists(a[["archivo"]])) { cat("Falta", a[["archivo"]], "- se omite", nm, "\n"); next }
  d <- as.data.frame(readRDS(file.path("data_clean", paste0(a[["datos"]], ".rds"))))
  pr <- leer_par(a[["archivo"]])
  r <- residuos_pit(d, pr$p, pr$z, T_exp, nsim = 2e5, semilla = 1)
  r$ajuste <- nm; r$z <- pr$z; res_pit[[nm]] <- r
  i <- which(d$dens == 10 & d$par == 2)
  cat(sprintf("\n== Residuos %s (z = %.2f): |r| > 3 en %d de %d; Shapiro-Wilk p = %.3f%s\n", nm, pr$z,
              sum(r$marca), nrow(r), stats::shapiro.test(r$residuo)$p.value,
              if (length(i)) sprintf("; ensayo 2/10: P(Y<=2) = %.1e, residuo = %.2f", r$p_baja[i], r$residuo[i]) else ""))
}
write.csv(do.call(rbind, res_pit), file.path("tablas", "diagnostico_residuos_pit.csv"), row.names = FALSE)
cat("\nArchivos: tablas/diagnostico_cribado_D2_D3.csv, tablas/diagnostico_cribado_ensayo1.csv,",
    "tablas/diagnostico_residuos_pit.csv\n")
