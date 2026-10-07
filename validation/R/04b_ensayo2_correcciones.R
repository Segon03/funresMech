# ============================================================
# R/04b_ensayo2_correcciones.R  (v2 - PARALELIZADO)
#
# v1 de este script corria los puntos pendientes uno por uno en un
# solo proceso (secuencial) - un error de Claude al escribirlo: a
# ~3h por punto con 1 solo core ocupado de 24 disponibles, terminar
# los hasta 17 puntos pendientes hubiera tardado 40-50h mas. Esta
# version manda TODAS las tareas pendientes juntas a future_lapply,
# igual que 02_ensayo1_simulacion.R / 04_ensayo2_datos_reales.R, asi
# se reparten entre los N workers disponibles y el tiempo total
# pasa a ser ~el tiempo de UNA tarea, no la suma de todas.
#
# Los 6 puntos que v1 ya alcanzo a calcular (checkpoints ya guardados
# en results/ensayo2/) NO se pierden ni se recalculan - el
# checkpointing por archivo sigue igual.
#
# Corrige los mismos 2 problemas que v1:
# 1) Los 2 puntos de perfil de D2 que no convergieron (z=2.30,
#    z=2.90 con NLL~868-871): se descartan (si siguen con ese valor
#    absurdo) y se recalculan.
# 2) El IC95% de D3 no cerraba por arriba (grilla original topaba
#    en z=3.0): se agregan candidatos hasta z=6.0, TODOS EN PARALELO
#    de una - no hace falta "parar apenas cruza el umbral" como en
#    v1, porque en paralelo calcular de mas no cuesta tiempo extra
#    (el limitante es la tarea mas lenta del lote, no la cantidad).
#
# Correr DESPUES de 04_ensayo2_datos_reales.R. Luego correr
# 04c_ensayo2_reensamblar.R.
# ============================================================

library(funresMech)
library(future)
library(future.apply)
source(file.path("fun", "betabinom_model.R"))
source(file.path("fun", "profile_utils.R"))

T_exp   <- 1
itermax <- 50
NP      <- 40
reltol  <- 1e-2
n_sim_profile <- 1500

lower_prof <- c(a = 0.001, h = 0.001, k = 0.5, s = 0.001)
upper_prof <- c(a = 2.0,   h = 0.5,   k = 5.0, s = 0.5)

data_D2 <- readRDS(file.path("data_clean", "D2.rds"))
data_D3 <- readRDS(file.path("data_clean", "D3.rds"))

tareas <- list()

# ---- D2: los 2 puntos que no convergieron ------------------------
for (z in c(2.30, 2.90)) {
  f <- file.path("results", "ensayo2", sprintf("perfil_D2_z%.2f.rds", z))
  if (file.exists(f)) {
    prev <- readRDS(f)$nll
    if (!is.na(prev) && prev > 300) {
      cat(sprintf("D2 z=%.2f: NLL previo=%.1f no convergido, se descarta\n", z, prev))
      file.remove(f)
    }
  }
  if (!file.exists(f)) {
    tareas[[length(tareas) + 1]] <- list(lab = "D2", z = z, data_spp = data_D2)
  }
}

# ---- D3: candidatos de extension hasta z=6.0, todos de una -------
for (z in seq(3.1, 6.0, by = 0.2)) {
  f <- file.path("results", "ensayo2", sprintf("perfil_D3_z%.2f.rds", z))
  if (!file.exists(f)) {
    tareas[[length(tareas) + 1]] <- list(lab = "D3", z = z, data_spp = data_D3)
  }
}

cat(sprintf("Tareas pendientes: %d\n", length(tareas)))

if (length(tareas) == 0) {
  cat("Nada pendiente - correr directamente R/04c_ensayo2_reensamblar.R\n")
} else {
  n_cores   <- parallel::detectCores()
  n_workers <- max(1, n_cores - 1)
  cat(sprintf("Paralelizando en %d workers (de %d cores detectados).\n", n_workers, n_cores))
  plan(multisession, workers = n_workers)

  correr_tarea <- function(t) {
    library(funresMech)
    out_file <- file.path("results", "ensayo2",
                           sprintf("perfil_%s_z%.2f.rds", t$lab, t$z))
    if (file.exists(out_file)) {
      return(sprintf("%s z=%.2f: ya existe, se salta", t$lab, t$z))
    }
    res <- DEoptim::DEoptim(
      fn = function(par) funresMech:::negloglik_fixed_z(par, t$z, t$data_spp, T_exp, n_sim_profile),
      lower = lower_prof, upper = upper_prof,
      control = DEoptim::DEoptim.control(itermax = itermax, NP = NP, reltol = reltol, trace = FALSE)
    )
    saveRDS(list(z = t$z, nll = res$optim$bestval), out_file)
    sprintf("%s z=%.2f: nll=%.3f", t$lab, t$z, res$optim$bestval)
  }

  mensajes <- future_lapply(tareas, correr_tarea, future.seed = TRUE)
  invisible(lapply(mensajes, cat, "\n"))
}

# ---- Chequeo final: ¿llego a cruzar el umbral para D3? ------------
archivos_D3 <- list.files(file.path("results", "ensayo2"),
                           pattern = "^perfil_D3_z.*\\.rds$", full.names = TRUE)
puntos_D3 <- lapply(archivos_D3, readRDS)
z_D3   <- vapply(puntos_D3, function(p) p$z, numeric(1))
nll_D3 <- vapply(puntos_D3, function(p) p$nll, numeric(1))
min_nll_D3 <- 132.21564831
threshold_D3 <- min_nll_D3 + qchisq(0.95, 1)
cruzo <- any(z_D3 > 2.9 & nll_D3 >= threshold_D3)

cat(sprintf("\nUmbral a cruzar (D3): %.4f\n", threshold_D3))
if (cruzo) {
  cat("El perfil de D3 ya cruzo el umbral en el rango extendido - el IC deberia cerrar.\n")
} else {
  cat("ADVERTENCIA: el perfil de D3 TODAVIA no cruza el umbral hasta z=6.0.\n",
      "Revisar tabla de z vs nll y decidir si extender mas o reportar el limite\n",
      "superior como no identificado dentro de un rango razonable.\n")
}

cat("\nListo. Ahora correr R/04c_ensayo2_reensamblar.R\n")
