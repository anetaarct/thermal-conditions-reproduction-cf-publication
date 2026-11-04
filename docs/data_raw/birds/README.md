# Collared flycatcher (CF) source datasets

## Source dataset used in the analysis lineage

`CF21980-2025.csv`: 15,984 records spanning seasons 1980–2025. Semicolon-delimited CSV, preserved byte-for-byte from the original GitHub upload.

- Source repository: https://github.com/anetaarct/gotland-flycatcher-phenology
- Source commit: `7c5344958c18a81e853838c0b7ad9fea751bec57` (28 January 2026).
- Original path at that commit: `CF21980-2025.csv`.
- SHA-256: `c7f9a9981cbcba9b065fa67e77e00ee40d5612a9fd990aa1fbcf1b2b157c62e6`.

The 14,693 records from 1980–2019 agree in LD, CS, HS, FLEDGE, NEWRECRUIT and FAGE_CAT with the sheet `Corrected for CS,FL,NEWR` in `CF_non-exp_breeding_events_1980-2019_RM_28.01.26.xlsx`. The file includes later-season records as well. Its import notes specify LD relative to 1 May, NEWRECRUIT rather than RECRUIT, and warn that HD contains errors. Missing observations were retained.

Data lineage: CF21980-2025.csv -> DB_MISMATCH.csv -> local DB_MISMATCH_current_onset_1980-2025.csv -> final trait-specific samples. Saved Git versions showed no FLEDGE changes across those transitions. NDVI and mismatch were not predictors in the final reproductive models.

For 2025, LD is observed in 198 source records and 179 records enter the final clutch-size sample. FLEDGE is not observed in 2025.

## Reproducing the models

Use `data_derived/analysis/clutch_model_data.csv` and `data_derived/analysis/fledging_model_data.csv` for the final reproductive fits. These include the final sample selection and prepared predictors; their README documents variables and centring. The source archives retain original ring identifiers and original field columns. They do not replace the final model inputs.
