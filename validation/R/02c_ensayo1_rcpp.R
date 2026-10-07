# ============================================================
# R/02c_ensayo1_rcpp.R
#
# ENSAYO 1 (definitivo): recuperacion de parametros conocidos con el
# MOTOR RCPP y la busqueda en ESCALA LOG (los mismos de 04h, Ensayo 2).
# Reemplaza al Ensayo 1 viejo (02/02b/03: motor R puro, n_sim 300,
# datos Beta-Binomial con densidades 5-100 y T = 1), que queda
# archivado como registro y NO va al manuscrito.
#
# Diseno (acordado el 02-oct-2026):
#   - Datos simulados con el PROPIO modelo de Okuyama (motor_okuyama_cpp,
#     modo "okuyama"), con el diseno experimental real de D3: mismas 13
#     densidades (2-50) y mismo numero de ensayos por densidad (64 en
#     total), T = 24 h.
#   - 3 escenarios de parametros verdaderos:
#       A  tipo II           z = 1     a = 0.05    h = 0.6    k = 1      s = 0.5
#       B  tipo III, "D3"    z = 2.5   a = 0.0107  h = 0.626  k = 0.821  s = 0.509
#       C  tipo III, "D2"    z = 2.75  a = 0.00238 h = 0.606  k = 0.0563 s = 0.641
#     (B y C: estimaciones de D3 y D2 con el motor Rcpp, 04h; en C z se
#      redondea de 2.8 a 2.75 para que caiga en la grilla.)
#     A: con a = 0.05 la proporcion parasitada baja de ~0.69 (N=2) a ~0.39
#     (N=50). La primera version (a = 0.65) saturaba las densidades bajas
#     (100% parasitado hasta N=8) y el piloto en la nube mostro que asi z
#     casi no se distingue por abajo (NLL en z=0.5 vs z=1: 0.24 de diferencia).
#     C prueba con datos simulados si, con k chico, z deja de estar
#     identificado por arriba (hipotesis surgida de D2 en el Ensayo 2).
#   - 20 replicas por escenario (60 datasets).
#   - Estimacion igual que en el Ensayo 2: perfil de z con DEoptim en
#     escala log (ajustar_z_fijo), z_hat = minimo del perfil, IC95% =
#     NLL <= NLL_min + qchisq(.95,1)/2. Grilla de z: 0.25 a 5.0, paso
#     0.25 (20 puntos).
#   - Configuracion REDUCIDA respecto del Ensayo 2 por costo:
#     n_sim = 5000, ITERMAX = 150 (NP = 40 y N_REEVAL = 5 iguales).
#     Ruido Monte Carlo de una evaluacion ~0.3 unidades de NLL (vs ~0.2
#     con 10000), chico frente al umbral de 1.92.
#
# Costo estimado: 1200 tareas (60 datasets x 20 valores de z), ~1 h
# cada una con 23 workers -> ~50 h de pared en la PC de computo.
#
# Checkpointing: cada dataset simulado y cada punto del perfil se guarda
# apenas termina. Si se corta, volver a correr retoma lo que falta. Los
# datasets se generan con semillas fijas (reproducibles) y se guardan
# ANTES de ajustar: re-correr nunca los cambia.
#
# REQUISITOS: los mismos de 04h (Rcpp + Rtools). Archivos necesarios:
#   fun/motor_okuyama.cpp, fun/motor_rcpp.R, data_clean/D3.rds
# Correr desde la carpeta del proyecto (.Rproj), idealmente como
# Background Job:
#   rstudioapi::jobRunScript("R/02c_ensayo1_rcpp.R", workingDir = getwd())
#
# Despues: R/03b_ensayo1_rcpp_resumen.R
#
# CORRECCION (04-oct-2026), despues de que la corrida se cortara a las
# ~35 h con "Error en load(file = index_file): archivo de entrada vacio":
#   - La causa era la cache de compilacion de Rcpp compartida entre
#     workers; se corrigio en fun/motor_rcpp.R (cache por proceso).
#   - Ademas, ahora un error en UNA tarea ya no cancela toda la corrida:
#     se reintenta una vez y, si vuelve a fallar, se registra y se sigue.
#     Al final se informa cuantas fallaron; volver a correr el script
#     las retoma.
#   - Cada resultado se guarda primero en un archivo temporal y despues
#     se renombra, asi un corte a mitad de escritura no deja un .rds roto.
#     Si al arrancar se encuentra un .rds ilegible, se renombra a
#     *.corrupto y esa tarea se vuelve a calcular.
# ============================================================

library(future)
library(future.apply)
source(file.path("fun", "motor_rcpp.R"))

# ---- Configuracion ----------------------------------------------------------
MODO_PRUEBA <- FALSE   # TRUE: corrida de humo de pocos minutos (1 replica por
                       # escenario, 3 valores de z, config minima; resultados en
                       # results/ensayo1_rcpp_PRUEBA/). Correrla antes de la larga.

ESCENARIOS <- list(
  A = c(a = 0.05,    h = 0.6,   z = 1,    k = 1,      s = 0.5),
  B = c(a = 0.0107,  h = 0.626, z = 2.5,  k = 0.821,  s = 0.509),
  C = c(a = 0.00238, h = 0.606, z = 2.75, k = 0.0563, s = 0.641)
)
N_REP    <- 20
T_exp    <- 24
N_SIM    <- 5000
ITERMAX  <- 150
NP       <- 40
N_REEVAL <- 5
Z_GRID   <- round(seq(0.25, 5.0, by = 0.25), 2)
SEMILLA_BASE <- 20261002
N_WORKERS <- max(1, parallel::detectCores() - 1)

if (MODO_PRUEBA) {
  N_REP <- 1; N_SIM <- 300; ITERMAX <- 3; N_REEVAL <- 2
  Z_GRID <- c(1, 2.5, 4)
}
dir_res <- file.path("results", if (MODO_PRUEBA) "ensayo1_rcpp_PRUEBA" else "ensayo1_rcpp")
dir.create(file.path(dir_res, "datos"), showWarnings = FALSE, recursive = TRUE)

# ---- Motor ------------------------------------------------------------------
cat("Compilando/cargando el motor C++ ...\n")
cargar_motor_rcpp()

# ---- Diseno experimental (el de D3) y datasets simulados ---------------------
f_d3 <- file.path("data_clean", "D3.rds")
if (!file.exists(f_d3)) stop("No se encontro ", f_d3, ". Correr 01_extract_forage_d2_d3.R primero.")
dens_diseno <- sort(as.data.frame(readRDS(f_d3))$dens)   # 64 ensayos, 13 densidades

simular_dataset <- function(p, dens, T, semilla) {
  set.seed(semilla)
  n_por_dens <- table(dens)
  do.call(rbind, lapply(names(n_por_dens), function(H) {
    H <- as.integer(H)
    y <- motor_okuyama_cpp(H, p[["a"]], p[["h"]], p[["z"]], p[["k"]], p[["s"]], T,
                           as.integer(n_por_dens[[as.character(H)]]), 0L, TRUE)
    data.frame(dens = H, par = y)
  }))
}

datasets <- list()
for (esc in names(ESCENARIOS)) for (r in seq_len(N_REP)) {
  id <- sprintf("%s_rep%02d", esc, r)
  f  <- file.path(dir_res, "datos", paste0(id, ".rds"))
  if (!file.exists(f)) {
    semilla <- SEMILLA_BASE + 1000 * match(esc, names(ESCENARIOS)) + r
    d <- simular_dataset(ESCENARIOS[[esc]], dens_diseno, T_exp, semilla)
    saveRDS(list(datos = d, escenario = esc, replica = r, verdad = ESCENARIOS[[esc]],
                 semilla = semilla, T = T_exp), f)
  }
  guardado <- readRDS(f)
  if (!isTRUE(all.equal(guardado$verdad, ESCENARIOS[[esc]])))
    stop("El dataset guardado ", f, " fue simulado con otros parametros verdaderos que los\n",
         "  de ESCENARIOS (", esc, "). Borrar ", file.path(dir_res, "datos", paste0(esc, "_rep*.rds")),
         " y las carpetas ", file.path(dir_res, paste0(esc, "_rep*")), " y volver a correr.")
  datasets[[id]] <- guardado$datos
}
cat(sprintf("Datasets: %d (%s x %d replicas), guardados en %s\n",
            length(datasets), paste(names(ESCENARIOS), collapse = "/"), N_REP,
            file.path(dir_res, "datos")))

# ---- Tareas -------------------------------------------------------------------
tareas <- list()
for (z in Z_GRID) for (id in names(datasets)) {   # z por fuera: los primeros
  dir.create(file.path(dir_res, id), showWarnings = FALSE)   # resultados ya cubren
  f <- file.path(dir_res, id, sprintf("perfil_z%.2f.rds", z)) # todos los datasets
  if (file.exists(f)) {
    ok <- tryCatch({ r <- readRDS(f); is.list(r) && is.finite(r$nll) }, error = function(e) FALSE)
    if (isTRUE(ok)) next
    cat("  Archivo ilegible, se recalcula:", f, "\n")
    file.rename(f, paste0(f, ".corrupto"))
  }
  tareas[[length(tareas) + 1]] <- list(id = id, z = z, archivo = f)
}
n_total <- length(datasets) * length(Z_GRID)
t_eval <- system.time(for (i in 1:3)
  nll_motor(ESCENARIOS$B[c("a", "h", "k", "s")], ESCENARIOS$B[["z"]],
            datasets[[1]], T_exp, N_SIM))[["elapsed"]] / 3
evals_por_tarea <- NP * (ITERMAX + 1) + N_REEVAL
n_workers <- max(1, min(N_WORKERS, length(tareas)))
horas <- length(tareas) * evals_por_tarea * t_eval / 3600
cat(sprintf(paste0(
  "Tareas: %d de %d pendientes (%d datasets x %d valores de z)\n",
  "1 evaluacion de la NLL (n_sim=%d): %.2f s (sin carga; con todos los workers\n",
  "  ocupados suele ser ~2.5x mas lenta) -> ~%.0f horas-nucleo\n",
  "Con %d workers: ~%.0f-%.0f horas de pared (estimacion)\n\n"),
  length(tareas), n_total, length(datasets), length(Z_GRID), N_SIM, t_eval, horas,
  n_workers, horas / n_workers, 2.5 * horas / n_workers))

if (length(tareas) == 0) {
  cat("No hay tareas pendientes. Correr R/03b_ensayo1_rcpp_resumen.R\n")
} else {
  plan(multisession, workers = n_workers)

  correr_tarea <- function(tarea) {
    if (file.exists(tarea$archivo)) return(sprintf("%s z=%.2f: ya existe", tarea$id, tarea$z))
    un_intento <- function() {
      source(file.path("fun", "motor_rcpp.R"))
      cargar_motor_rcpp()
      res <- ajustar_z_fijo(tarea$z, datasets[[tarea$id]], T = CFG$T_exp, nsim = CFG$N_SIM,
                            escala = "log", itermax = CFG$ITERMAX, NP = CFG$NP,
                            n_reeval = CFG$N_REEVAL)
      res$dataset <- tarea$id
      tmp <- paste0(tarea$archivo, ".tmp", Sys.getpid())   # escritura atomica
      saveRDS(res, tmp)
      if (!file.rename(tmp, tarea$archivo)) stop("no se pudo renombrar ", tmp)
      res
    }
    res <- tryCatch(un_intento(), error = function(e1) {
      Sys.sleep(5)
      tryCatch(un_intento(), error = function(e2) e2)
    })
    if (inherits(res, "error"))
      return(sprintf("ERROR %s z=%.2f: %s", tarea$id, tarea$z, conditionMessage(res)))
    sprintf("%s z=%.2f: nll=%.3f (a=%.3g h=%.3f k=%.3g s=%.3f) %.1f min%s",
            tarea$id, tarea$z, res$nll, res$par[["a"]], res$par[["h"]], res$par[["k"]],
            res$par[["s"]], res$segundos / 60,
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
  n_err <- sum(startsWith(unlist(mensajes), "ERROR"))
  if (n_err > 0) cat(sprintf(paste0("\nATENCION: %d tareas fallaron (ver lineas 'ERROR' arriba). ",
                                    "Volver a correr este script las retoma.\n"), n_err))
  cat(sprintf("\nListo en %.1f horas. Correr R/03b_ensayo1_rcpp_resumen.R%s\n",
              as.numeric(difftime(Sys.time(), t0, units = "hours")),
              if (MODO_PRUEBA) " (con MODO_PRUEBA <- TRUE tambien alli)" else ""))
  plan(sequential)
}
