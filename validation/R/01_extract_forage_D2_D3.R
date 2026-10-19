# ============================================================================
# 01_extract_forage_D2_D3.R
#
# Extracts datasets D2 (FoRAGE curve 57) and D3 (curve 58) from the FoRAGE
# database, version 5 (Uiterwaal et al. 2022, Ecology 103: e3706). They are
# the second- and third-instar panels of the first figure of Kishani Farahani
# and Goldansaz (2013, Eur J Entomol 110: 277-283): Apanteles myeloenta
# attacking Ectomyelois ceratoniae for 24 h. The curves were digitized point
# by point from the published figure, so trials with the same count overlap
# and mostly the distinct values at each density are present (S1 Text,
# section 11).
#
# Requires validation/data_raw/FoRAGE_db_V5_Dec_20_2024_original_curves.csv,
# which is not redistributed; obtain it from the FoRAGE data paper.
# Output: validation/data/D2.csv, validation/data/D3.csv (columns dens, par),
# and validation/data/forage_D2_D3_raw.csv (with metadata).
# ============================================================================

csv_path <- file.path("validation", "data_raw", "FoRAGE_db_V5_Dec_20_2024_original_curves.csv")
if (!file.exists(csv_path))
  stop("File not found: ", csv_path, ". Download the FoRAGE database (version 5, original curves) ",
       "and place it in validation/data_raw/.")

curves <- utils::read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
names(curves)[1] <- "data_set_id"        # the first column (curve identifier) has no name in the file

ids <- c(D2 = 57, D3 = 58)
raw <- curves[curves$data_set_id %in% ids, ]
raw <- data.frame(dataset = names(ids)[match(raw$data_set_id, ids)], data_set_id = raw$data_set_id,
                  source = raw$Source, predator = raw$Predator, prey = raw$Prey,
                  dens = round(raw[["Original x"]]), par = round(raw[["Original y"]]),
                  T_hours = raw[["Trial duration (h)"]],
                  sample_size = raw[["Sample size (per density)"]])

stopifnot("counts larger than the host density" = all(raw$par <= raw$dens),
          "trial duration different from 24 h" = all(raw$T_hours == 24),
          "no rows found for D2 and D3" = nrow(raw) > 0)

dir.create(file.path("validation", "data"), showWarnings = FALSE)
for (nm in names(ids)) {
  d <- raw[raw$dataset == nm, c("dens", "par")]
  utils::write.csv(d, file.path("validation", "data", paste0(nm, ".csv")), row.names = FALSE)
  cat(sprintf("%s: %d trials written to validation/data/%s.csv\n", nm, nrow(d), nm))
}
utils::write.csv(raw, file.path("validation", "data", "forage_D2_D3_raw.csv"), row.names = FALSE)
