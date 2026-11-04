# Original collared flycatcher (CF) source dataset

`CF_adult_birds_1980-2019.csv` contains 23,911 original adult-bird breeding records from the Gotland, Sweden, nest-box monitoring programme, covering 1980–2019. The CSV was copied unchanged from the project source archive.

SHA-256: `e1d7e3a6827b8c3d51afa26f6a8a4e3c06086c9ab7f3ba321fb53f583a5f6612`

The checked fields (ID, YEAR, FRING, LD, CS, HS and FLEDGE) agree across all records with the original workbook `RM_CF_adults_1980-2019_14.07.23.xlsx`, after normalising numeric formatting and missing values. Ring identifiers and the original field columns are retained in this source archive. The check establishes agreement in these fields, rather than bytewise identity of the Excel and CSV files.

## Relationship to the publication analyses

This original file is not the final model input. The reproductive preparation workflow used the later `DB_MISMATCH_current_onset_1980-2025.csv` dataset. That version contains 15,984 records, includes post-2019 seasons and differs from this archive. Matching by YEAR, FRING, AREA and BOX identified 14,692 unambiguous record matches and 312 FLEDGE differences: 247 missing-to-observed, 57 observed-to-different-observed and eight observed-to-missing. Most occurred in 2018–2019. A previous ID/YEAR-only comparison was unreliable because some IDs were missing. A complete change history is not available in the published materials.

Use `data_derived/analysis/clutch_model_data.csv` and `data_derived/analysis/fledging_model_data.csv` to reproduce the final models. These contain the final sample selection and prepared predictors. Their README documents variables and centring. The final models did not use NDVI or phenological mismatch as predictors.

## Selected original fields

- ID: original record identifier; YEAR: breeding season.
- FRING: original female ring identifier; MRING: original male ring identifier.
- LD: laying date in the source coding; CS: clutch size; HS: hatchling number; FLEDGE: fledgling number.
- Other field columns are preserved as supplied. Their original coding should not be reinterpreted without the field-data documentation.

Data provenance and the distinction between the original archive, the later breeding dataset and the final model inputs are described on the website's A1 page.
