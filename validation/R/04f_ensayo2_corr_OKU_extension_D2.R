# ============================================================
# R/04f_ensayo2_corr_OKU_extension_D2.R
#
# Extiende el perfil de z de D2 (corr_OKU) desde z = 3.1 hasta
# z = 5.9 (paso 0.2, 15 puntos nuevos), con EXACTAMENTE la misma
# configuracion que 04d_ensayo2_corr_OKU.R.
#
# POR QUE: en 04d la grilla de D2 quedo en seq(0.5, 3.0, 0.2), es
# decir, termina en z = 2.9. Pero:
#   - Okuyama (2026) evalua z de 0.1 a 5 (Sec. 2.3) y reporta para
#     D2 (simulacion) z_hat = 3.1 e IC = [1.1, 4.6] (Tabla S1): tanto
#     su estimador como su limite superior quedan FUERA de nuestra
#     grilla.
#   - Con la grilla actual, el perfil de D2 nunca cruza el umbral
#     del lado derecho -> el limite superior del IC sale NA
#     (censurado: solo sabemos que es >= 2.9).
# Se extiende hasta 5.9 para que D2 y D3 queden con la misma grilla
# (0.5 a 5.9, paso 0.2).
#
# NOVEDAD RESPECTO DE 04d: cada punto del perfil guarda tambien los
# parametros (a, h, k, s) del optimo condicional, no solo z y NLL.
# Asi, si el minimo del perfil de D2 cae en uno de estos puntos
# nuevos, se pueden reportar a, h, k, s en ese minimo (como hace
# Okuyama) y re-evaluar la NLL con n_sim alto si hiciera falta.
#
# COSTO: 15 tareas -> una sola "tanda" con 15 workers. Tiempo de
# pared esperado ~1-1.5 dias (una tanda de 04d tardo ~30-36 h).
#
# NO borra ni modifica ningun archivo de 04d: solo agrega archivos
# perfil_D2_corrOKU_z3.10.rds ... z5.90.rds en la misma carpeta.
# Checkpointing igual que 04d: si se corta, correr de nuevo retoma
# solo lo que falta.
#
# IMPORTANTE: correr desde la carpeta del proyecto (abrir el .Rproj).
# ============================================================

library(funresMech)
library(future)
library(future.apply)

res_dir <- file.path("results", "ensayo2_corr_OKU")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Configuracion: IDENTICA a 04d_ensayo2_corr_OKU.R -------------------
# (no cambiar ninguno de estos valores: los puntos nuevos tienen que ser
# comparables con los 13 puntos de D2 que ya existen)
T_exp   <- 24
itermax <- 200
NP      <- 40
N_SIM   <- 3000
LIMITE_SEG_POR_EVAL <- 180
lower_ok <- c(a = 1e-6, h = 0.01, k = 1e-6, s = 0)
upper_ok <- c(a = 5,    h = 5,    k = 5,    s = 3)

evaluar_con_limite <- function(expr_fn, limite_seg = LIMITE_SEG_POR_EVAL) {
  setTimeLimit(elapsed = limite_seg, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE))
  tryCatch(expr_fn(), error = function(e) 1e10)
}

# ---- Puntos nuevos --------------------------------------------------------
z_nuevos <- round(seq(3.1, 5.9, by = 0.2), 2)

data_path <- file.path("data_clean", "D2.rds")
if (!file.exists(data_path)) stop("No se encontro ", data_path, ". Correr 01_extract_forage_d2_d3.R primero.")
data_D2 <- readRDS(data_path)
stopifnot(nrow(data_D2) == 73)   # mismo n que el Script S3 de Okuyama

pendientes <- z_nuevos[!file.exists(file.path(res_dir, sprintf("perfil_D2_corrOKU_z%.2f.rds", z_nuevos)))]
cat(sprintf("Extension del perfil de D2: %d puntos en total, %d pendientes.\n",
            length(z_nuevos), length(pendientes)))
if (length(pendientes) == 0) {
  cat("Nada que correr. Pasar a 04g_ensayo2_corr_OKU_reensamblar_v2.R.\n")
} else {
  n_cores   <- parallel::detectCores()
  n_workers <- max(1, min(length(pendientes), n_cores - 1))
  cat(sprintf("Paralelizando en %d workers (de %d cores detectados).\n", n_workers, n_cores))
  cat("OJO: future_lapply no imprime nada hasta que terminan TODAS las tareas.\n",
      "Para ver el avance, mirar como aparecen los archivos perfil_D2_corrOKU_z*.rds\n",
      "en results/ensayo2_corr_OKU/.\n", sep = "")
  plan(multisession, workers = n_workers)

  correr_punto <- function(z) {
    library(funresMech)
    out_file <- file.path(res_dir, sprintf("perfil_D2_corrOKU_z%.2f.rds", z))
    if (file.exists(out_file)) return(sprintf("D2 z=%.2f: ya existe, se salta", z))
    t0 <- Sys.time()
    res <- DEoptim::DEoptim(
      fn = function(par) evaluar_con_limite(function() {
        funresMech:::negloglik_fixed_z(par, z_fixed = z, data_spp = data_D2,
                                        T = T_exp, n_sim = N_SIM)
      }),
      lower = lower_ok, upper = upper_ok,
      control = DEoptim::DEoptim.control(itermax = itermax, NP = NP, trace = FALSE)
    )
    par <- res$optim$bestmem
    names(par) <- c("a", "h", "k", "s")
    saveRDS(list(z = z, nll = res$optim$bestval, par = par,
                 n_sim = N_SIM, itermax = itermax, NP = NP, T = T_exp,
                 horas = as.numeric(Sys.time() - t0, units = "hours")),
            out_file)
    sprintf("D2 z=%.2f: nll=%.3f (%.1f h)", z, res$optim$bestval,
            as.numeric(Sys.time() - t0, units = "hours"))
  }

  mensajes <- future_lapply(pendientes, correr_punto, future.seed = TRUE)
  invisible(lapply(mensajes, cat, "\n"))
  plan(sequential)
}

cat("\nListo. Correr 04g_ensayo2_corr_OKU_reensamblar_v2.R.\n")
