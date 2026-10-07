# ============================================================
# R/04j_ensayo2_rcpp_extension_D2.R
#
# Extension del perfil de z de D2 (escala "log") MAS ALLA de la
# grilla de Okuyama (0.5-6.0), para ver hasta donde llega el limite
# superior del IC95%.
#
# Por que (resultado de 04h/04i, corrida del 30-sep-2026):
#   - D2, escala log: entre z ~1.6 y 6.0 el perfil queda < 0.6
#     unidades de NLL por encima del minimo (umbral IC95 = +1.92),
#     sin ningun parametro en la cota. El IC superior queda
#     censurado en la grilla: [1.2, >=6.0] (Okuyama: [1.1, 4.6]).
#     A medida que z sube, 'a' baja casi exactamente como H_ref^-z
#     (compensacion a-z) y el ajuste casi no empeora.
#   - D2, escala okuyama (busqueda lineal de a en [0,5], como el
#     Script S3): desde z = 2.5 'a' queda en la cota inferior
#     (a < 0.005) y desde z ~4.3 la NLL salta de forma erratica y la
#     NLL re-evaluada explota (705, 564, ...). Ese deterioro de la
#     busqueda es lo que hace cruzar el umbral cerca de 4.2-4.6:
#     el limite superior de Okuyama (4.6) seria un artefacto de la
#     busqueda, no informacion de los datos.
#
# Esta extension responde la pregunta que queda: con la busqueda
# log (que no tiene ese problema: lambda_ref = a*H_ref^z tiene cotas
# fijas que no dependen de z), el perfil de D2
#   (a) se aplana del todo  -> z no identificable por arriba,
#                              IC informado como [1.2, Inf), o
#   (b) cruza +1.92 en algun z alto -> ese es el limite superior.
#
# Misma configuracion que 04h (N_SIM 10000, ITERMAX 200, NP 40,
# N_REEVAL 5, motor Rcpp, corte por saturacion). Solo cambia la
# grilla de z. Los resultados van a una carpeta APARTE
# (results/ensayo2_rcpp/log_ext/) para no mezclar la extension con
# la grilla 0.5-6.0 que se compara 1 a 1 con Okuyama; 04i los lee
# y los analiza en una seccion propia.
#
# Costo: 7 tareas en paralelo (una por valor de z). Con los tiempos
# de 04h (~2-2.5 h por tarea con la PC cargada) deberia tardar
# ~2-2.5 h de pared. Checkpointing igual que 04h: si se corta, volver
# a correr retoma lo que falta.
#
# Correr desde la carpeta del proyecto (.Rproj). Recomendado como
# Background Job:  rstudioapi::jobRunScript("R/04j_ensayo2_rcpp_extension_D2.R")
#
# Despues: R/04i_ensayo2_rcpp_reensamblar.R (detecta log_ext/ solo)
# ============================================================

library(future)
library(future.apply)
source(file.path("fun", "motor_rcpp.R"))

# ---- Configuracion ----------------------------------------------------------
MODO_PRUEBA <- FALSE   # TRUE: corrida de humo de pocos minutos (config reducida,
                       # resultados en results/ensayo2_rcpp_PRUEBA/log_ext/)

ESCALA   <- "log"      # solo la principal: en escala "okuyama" la busqueda ya
                       # se degrada desde z ~4.3, extenderla no informa nada
LAB      <- "D2"
T_exp    <- 24
N_SIM    <- 10000      # = 04h
ITERMAX  <- 200        # = 04h
NP       <- 40         # = 04h
N_REEVAL <- 5          # = 04h
Z_EXT    <- c(6.5, 7, 8, 10, 12, 15, 20)
N_WORKERS <- max(1, parallel::detectCores() - 1)

if (MODO_PRUEBA) {
  N_SIM <- 500; ITERMAX <- 3; N_REEVAL <- 2
  Z_EXT <- c(8, 20)
}
dir_res <- file.path("results", if (MODO_PRUEBA) "ensayo2_rcpp_PRUEBA" else "ensayo2_rcpp", "log_ext")
dir.create(dir_res, showWarnings = FALSE, recursive = TRUE)

# ---- Motor y datos ------------------------------------------------------------
cat("Compilando/cargando el motor C++ ...\n")
cargar_motor_rcpp()

f_dat <- file.path("data_clean", paste0(LAB, ".rds"))
if (!file.exists(f_dat)) stop("No se encontro ", f_dat, ". Correr 01_extract_forage_d2_d3.R primero.")
datasets <- setNames(list(as.data.frame(readRDS(f_dat))), LAB)

# Chequeo numerico rapido a z alto (que la NLL sea finita con parametros
# tipicos del extremo de la grilla de 04h: en z = 6.0, lambda_ref =
# a*H_ref^z = 10.8, h = 0.63, k = 0.017, s = 0.66).
# Nota: lambda_ref crece con z (5.0 en z = 2.8, 10.8 en z = 6.0); la cota
# superior de la busqueda log es 1e3. Si a z = 15-20 llegara a la cota,
# ajustar_z_fijo() lo marca (en_cota) y 04i lo informa.
H_ref <- h_ref_de(datasets[[LAB]])
par_prueba <- c(a = 10.8 / H_ref^max(Z_EXT), h = 0.63, k = 0.017, s = 0.66)
nll_prueba <- nll_motor(par_prueba, max(Z_EXT), datasets[[LAB]], T_exp, 500)
cat(sprintf("Chequeo a z = %g: NLL = %.2f (%s)\n", max(Z_EXT), nll_prueba,
            if (is.finite(nll_prueba)) "ok" else "NO FINITA - revisar antes de seguir"))
if (!is.finite(nll_prueba)) stop("La NLL no es finita a z alto.")

# ---- Tareas -------------------------------------------------------------------
tareas <- list()
for (z in Z_EXT) {
  f <- file.path(dir_res, sprintf("perfil_%s_z%.2f.rds", LAB, z))
  if (!file.exists(f)) tareas[[length(tareas) + 1]] <- list(escala = ESCALA, lab = LAB, z = z, archivo = f)
}
t_eval <- system.time(for (i in 1:3)
  nll_motor(c(0.02, 0.5, 1, 0.5), 2, datasets[[LAB]], T_exp, N_SIM))[["elapsed"]] / 3
evals_por_tarea <- NP * (ITERMAX + 1) + N_REEVAL
n_workers <- max(1, min(N_WORKERS, length(tareas)))
cat(sprintf(paste0(
  "\nTareas: %d de %d pendientes (%s, %s, z = %s)\n",
  "1 evaluacion de la NLL (n_sim=%d): %.2f s -> ~%.1f h por tarea (estimacion;\n",
  "  a z alto puede variar). Con %d workers en paralelo: ~%.1f h de pared.\n\n"),
  length(tareas), length(Z_EXT), ESCALA, LAB, paste(Z_EXT, collapse = ", "),
  N_SIM, t_eval, evals_por_tarea * t_eval / 3600, n_workers,
  ceiling(length(tareas) / n_workers) * evals_por_tarea * t_eval / 3600))

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
    sprintf("%-4s %s z=%5.2f: nll=%.3f reeval=%.3f (a=%.3g h=%.3f k=%.3g s=%.3f) %.1f min%s",
            tarea$escala, tarea$lab, tarea$z, res$nll, res$nll_reeval,
            res$par[["a"]], res$par[["h"]], res$par[["k"]], res$par[["s"]], res$segundos / 60,
            if (any(res$en_cota)) paste0("  [en cota: ", paste(names(res$en_cota)[res$en_cota], collapse = ","), "]") else "")
  }

  CFG <- list(T_exp = T_exp, N_SIM = N_SIM, ITERMAX = ITERMAX, NP = NP, N_REEVAL = N_REEVAL)
  t0 <- Sys.time()
  mensajes <- future_lapply(
    tareas, correr_tarea,
    future.seed       = TRUE,
    future.scheduling = Inf,
    future.globals    = list(datasets = datasets, CFG = CFG),
    future.packages   = c("Rcpp", "DEoptim")
  )
  invisible(lapply(mensajes, function(m) cat(m, "\n")))
  cat(sprintf("\nListo en %.1f horas. Correr R/04i_ensayo2_rcpp_reensamblar.R%s\n",
              as.numeric(difftime(Sys.time(), t0, units = "hours")),
              if (MODO_PRUEBA) " (con MODO_PRUEBA <- TRUE tambien alli)" else ""))
  plan(sequential)
}
