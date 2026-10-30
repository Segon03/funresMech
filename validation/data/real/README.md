# Published datasets used in 25_real_data.R (besides FoRAGE)

- `DAuro2020_RAW_DATA.csv`: raw data of D'Auro F, Luna MG, Liljesthröm GG (2020). Functional
  response and egg distribution among hosts by *Pseudapanteles dignus* (Hymenoptera: Braconidae),
  a larval endoparasitoid of *Tuta absoluta* (Lepidoptera: Gelechiidae). Zenodo, record 4105529,
  CC BY 4.0. One row per host: `treatment` (hosts offered to one female for 24 h), `replicate`,
  `individual` (host within the replicate), `value` (parasitoid eggs found by dissection).
- `DAuro2020_trials.csv`: the same data by trial (built from the raw file): hosts offered (`dens`),
  hosts recovered and dissected, hosts with at least one egg (`par`) and eggs. Some hosts were not
  recovered (54 of 836); `par` counts parasitized hosts among those recovered.
- `NunezCampero2025.csv`: Núñez-Campero et al. (2025), Insects 16: 919. Four parasitoids of
  *Anastrepha fraterculus*; one female per trial, 24 h, hosts not replaced; parasitism scored by
  adult emergence (`par`), with dead pupae in `dead_flies`.
- The *Gammarus* data (Paterson et al. 2015) are read from the frair package (`data(gammarus)`).

The FoRAGE database (DeLong and Uiterwaal, KNB, doi:10.5063/F1RX99KB, version 5, December 2024)
is not redistributed: place `FoRAGE_db_V5_Dec_20_2024_original_curves.csv` and
`FoRAGE_db_V5_Dec_20_2024_sources_and_meta.csv` in `validation/data_raw/`.
