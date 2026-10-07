# ============================================================
# R/04h_ensayo2_rcpp_log.R
#
# Ensayo 2 (datos reales D2 y D3) con el MOTOR RCPP y la
# REPARAMETRIZACION LOG. Reemplaza a 04d/04f (motor R puro de
# funresMech v1.0.4) como corrida de validacion definitiva.
#
# Que cambia respecto de 04d/04f (detalle y justificacion en
# docs/MOTOR_RCPP_cambios_y_validacion.md):
#
#  1. Motor de simulacion en C++ (fun/motor_okuyama.cpp), el mismo
#     proceso que el Script S5 de Okuyama (verificado bit a bit en
#     R/00b_motor_rcpp_pruebas.R). ~65-80x mas rapido que el R puro.
#  2. Parametrizacion de s como Okuyama: s = desvio estandar NATURAL
#     del tiempo de manipulacion, buscado como s = u^2 (s.sqrt = TRUE
#     del Script S3). Asi a, h, k, s son directamente comparables con
#     la Tabla S1.
#  3. N_SIM = 10000, igual que Okuyama (04d/04f usaban 3000 por costo).
#  4. Corte por saturacion (exacto): elimina la trampa h~0/k~0 sin
#     poner pisos artificiales a h ni a k (04d/04f usaban h >= 0.01).
#  5. Busqueda de DEoptim en escala log (escala = "log"):
#     log(a*H_ref^z), log h, log k, sqrt(s). Resuelve la mala
#     exploracion de a y k chicos a z alto que dejaba el perfil
#     "serruchado" en 04f.
#  6. Grilla de z con paso 0.1 (como Okuyama), de 0.5 a 6.0.
#  7. Se guarda en cada punto: parametros, NLL de DEoptim (bestval,
#     comparable con Okuyama) y NLL re-evaluada con semillas nuevas
#     (sin el sesgo optimista del minimo de evaluaciones ruidosas).
#
# Opcional (ESCALAS de abajo): correr ademas la busqueda EXACTA de
# Okuyama (escala = "okuyama": lineal, [0,5]^4) con el mismo motor.
# Sirve como "ablacion" para el paper: muestra cuanto aporta la
# reparametrizacion log. Las tareas "log" se despachan primero.
#
# REQUISITOS: Rcpp y Rtools (Windows) - ver docs/. Correr ANTES
# R/00b_motor_rcpp_pruebas.R y verificar que todas las pruebas pasen.
#
# Checkpointing: cada punto del perfil se guarda en su propio .rds
# apenas termina; si se corta, volver a correr retoma lo que falta.
# Correr desde la carpeta del proyecto (.Rproj).
#
# Despues: R/04i_ensayo2_rcpp_reensamblar.R
# ============================================================

library(future)
library(future.apply)
source(file.path("fun", "motor_rcpp.R"))

# ---- Configuracion ----------------------------------------------------------
MODO_PRUEBA <- FALSE   # TRUE: corrida de humo de pocos minutos (config reducida,
                       # resultados en results/ensayo2_rcpp_PRUEBA/) para verificar
                       # que todo funciona antes de lanzar la corrida larga.

ESCALAS  <- c("log", "okuyama")   # "log" = principal; "okuyama" = ablacion
                                  # (sacar "okuyama" si se quiere solo la principal:
                                  #  reduce el tiempo a la mitad)
T_exp    <- 24
N_SIM    <- 10000      # Okuyama (2026), Script S3
ITERMAX  <- 200        # Okuyama (2026), Script S3
NP       <- 40         # default de DEoptim para 4 parametros (10 x 4), como Okuyama
N_REEVAL <- 5          # re-evaluaciones de la NLL en el optimo (semillas nuevas)
Z_GRID   <- round(seq(0.5, 6.0, by = 0.1), 2)   # paso 0.1 como Okuyama
labs     <- c("D2", "D3")
N_WORKERS <- max(1, parallel::detectCores() - 1)

if (MODO_PRUEBA) {
  N_SIM <- 500; ITERMAX <- 3; N_REEVAL <- 2
  Z_GRID <- c(1.5, 2.5, 3.5)
}
dir_res <- file.path("results", if (MODO_PRUEBA) "ensayo2_rcpp_PRUEBA" else "ensayo2_rcpp")
for (esc in ESCALAS) dir.create(file.path(dir_res, esc), showWarnings = FALSE, recursive = TRUE)

# ---- Motor y datos ------------------------------------------------------------
cat("Compilando/cargando el motor C++ ...\n")
cargar_motor_rcpp()   # compila UNA vez aca; los workers reutilizan la cache

datasets <- setNames(lapply(labs, function(lab) {
  f <- file.path("data_clean", paste0(lab, ".rds"))
  if (!file.exists(f)) stop("No se encontro ", f, ". Correr 01_extract_forage_d2_d3.R primero.")
  as.data.frame(readRDS(f))
}), labs)

# ---- Benchmark ----------------------------------------------------------------
t_eval <- system.time(for (i in 1:3)
  nll_motor(c(0.02, 0.5, 1, 0.5), 2, datasets[["D2"]], T_exp, N_SIM))[["elapsed"]] / 3

# ---- Tareas ---------------------------------------------------------------------
tareas <- list()
for (esc in ESCALAS) for (lab in labs) for (z in Z_GRID) {
  f <- file.path(dir_res, esc, sprintf("perfil_%s_z%.2f.rds", lab, z))
  if (!file.exists(f)) tareas[[length(tareas) + 1]] <- list(escala = esc, lab = lab, z = z, archivo = f)
}
n_total <- length(ESCALAS) * length(labs) * length(Z_GRID)
evals_por_tarea <- NP * (ITERMAX + 1) + N_REEVAL
n_workers <- max(1, min(N_WORKERS, length(tareas)))
horas <- length(tareas) * evals_por_tarea * t_eval / 3600

cat(sprintf(paste0(
  "\n1 evaluacion de la NLL (D2, n_sim=%d): %.2f s\n",
  "Tareas: %d de %d pendientes (%s x %s x %d valores de z)\n",
  "~%d evaluaciones por tarea -> ~%.1f horas-nucleo en total\n",
  "Con %d workers: ~%.1f horas de pared (estimacion; D3 y z altos pueden variar)\n\n"),
  N_SIM, t_eval, length(tareas), n_total, paste(ESCALAS, collapse = "+"),
  paste(labs, collapse = "+"), length(Z_GRID), evals_por_tarea, horas,
  n_workers, horas / n_workers))

if (length(tareas) == 0) {
  cat("No hay tareas pendientes. Correr R/04i_ensayo2_rcpp_reensamblar.R\n")
} else {
  plan(multisession, workers = n_workers)

  correr_tarea <- function(tarea) {
    if (file.exists(tarea$archivo)) return(sprintf("%s %s z=%.2f: ya existe", tarea$escala, tarea$lab, tarea$z))
    source(file.path("fun", "motor_rcpp.R"))
    cargar_motor_rcpp()
    res <- ajustar_z_fijo(tarea$z, datasets[[tarea$lab]], T = CFG$T_exp, nsim = CFG$N_SIM,
                          escala = tarea$escala, itermax = CFG$ITERMAX, NP = CFG$NP,
                          n_reeval = CFG$N_REEVAL)
    res$dataset <- tarea$lab
    saveRDS(res, tarea$archivo)
    sprintf("%-7s %s z=%.2f: nll=%.3f reeval=%.3f (a=%.3g h=%.3f k=%.3g s=%.3f) %.1f min%s",
            tarea$escala, tarea$lab, tarea$z, res$nll, res$nll_reeval,
            res$par[["a"]], res$par[["h"]], res$par[["k"]], res$par[["s"]], res$segundos / 60,
            if (any(res$en_cota)) paste0("  [en cota: ", paste(names(res$en_cota)[res$en_cota], collapse = ","), "]") else "")
  }

  CFG <- list(T_exp = T_exp, N_SIM = N_SIM, ITERMAX = ITERMAX, NP = NP, N_REEVAL = N_REEVAL)
  t0 <- Sys.time()
  mensajes <- future_lapply(
    tareas, correr_tarea,
    future.seed       = TRUE,          # semillas L'Ecuyer independientes y reproducibles
    future.scheduling = Inf,           # una tarea por vez y en orden: "log" sale primero
    future.globals    = list(datasets = datasets, CFG = CFG),
    future.packages   = c("Rcpp", "DEoptim")
  )
  invisible(lapply(mensajes, function(m) cat(m, "\n")))
  cat(sprintf("\nListo en %.1f horas. Correr R/04i_ensayo2_rcpp_reensamblar.R%s\n",
              as.numeric(difftime(Sys.time(), t0, units = "hours")),
              if (MODO_PRUEBA) " (con MODO_PRUEBA <- TRUE tambien alli)" else ""))
  plan(sequential)
}
