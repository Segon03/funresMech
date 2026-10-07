# ============================================================
# R/01_extract_forage_d2_d3.R  (Seccion 4.1 del plan)
#
# Extrae D2 (FoRAGE 'Data set' = 57) y D3 ('Data set' = 58) desde
# el CSV crudo de FoRAGE v5, y los deja listos en data_clean/ con
# las columnas EXACTAS que espera funresMech:::fit_full(): dens, par.
#
# Requiere: data_raw/FoRAGE_db_V5_Dec_20_2024_original_curves.csv
# (ya copiado desde la carpeta Analisis_FoRAGE al crear el proyecto).
# ============================================================

library(readr)
library(dplyr)

csv_path <- file.path("data_raw", "FoRAGE_db_V5_Dec_20_2024_original_curves.csv")
if (!file.exists(csv_path)) {
  stop("No se encontro ", csv_path, ". Copiar el CSV de FoRAGE v5 ",
       "(original_curves) a data_raw/ antes de correr este script.")
}

curves <- read_csv(csv_path, show_col_types = FALSE)

# La primera columna del CSV (el ID de "Data set" que usan Okuyama 2026 y
# el data-paper de Uiterwaal et al. 2022 para identificar cada curva) NO
# tiene nombre en el encabezado del CSV original - readr la nombra
# automaticamente "...1". Se renombra aca para que el resto del script
# sea legible y no dependa de ese nombre autogenerado.
names(curves)[1] <- "data_set_id"

d_ids <- c(D2 = 57, D3 = 58)

forage_raw <- curves %>%
  filter(data_set_id %in% d_ids) %>%
  mutate(
    dataset_label = names(d_ids)[match(data_set_id, d_ids)],
    dens = round(`Original x`),   # columna que espera funresMech:::fit_full()
    par  = round(`Original y`)    # idem
  ) %>%
  select(dataset_label, data_set_id, Source, Predator, Prey,
         dens, par, T_hours = `Trial duration (h)`,
         sample_size = `Sample size (per density)`)

# Chequeos de sanidad ---------------------------------------------------
stopifnot(
  "hay conteos de parasitismo mayores a la densidad ofrecida" =
    all(forage_raw$par <= forage_raw$dens),
  "la duracion del ensayo no es de 24 h en todas las filas" =
    all(forage_raw$T_hours == 24),
  "no se encontraron filas para D2 y D3" =
    nrow(forage_raw) > 0
)

dir.create("data_clean", showWarnings = FALSE)

for (lab in names(d_ids)) {
  df <- forage_raw %>% filter(dataset_label == lab) %>% select(dens, par)
  saveRDS(df, file.path("data_clean", paste0(lab, ".rds")))
  cat(sprintf("%s: %d filas guardadas en data_clean/%s.rds\n", lab, nrow(df), lab))
}

write_csv(forage_raw, file.path("data_clean", "forage_D2_D3_raw.csv"))
cat("\nListo. Ver tambien data_clean/forage_D2_D3_raw.csv (con metadatos, para trazabilidad).\n")
