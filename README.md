# Warmer springs and reproductive performance in collared flycatchers

Publication materials for *Larger clutches but lower egg-to-fledging success in warmer springs in a migratory passerine*.

**Results website:** https://anetaarct.github.io/thermal-conditions-reproduction-cf-publication/

This repository starts with a clean publication history. It contains the final maximum-likelihood analyses reported in the manuscript, their exported results and the companion Quarto website. Working drafts and superseded analyses are kept separately in a private archive.

## Reported models

Two distinct model sets answer separate questions. Temperature, laying-date and age associations retain the original temperature-association models below. Temporal trends use refits without temperature in either location or scale, retaining age, relative laying date, the same records, distributions and random-effect structures. Temporal estimates are not adjusted for temperature but remain standardised for age and relative laying date. They do not identify a causal warming effect.


- **Clutch size:** mean-parameterised COM-Poisson mixed model fitted by ML in glmmTMB. In the temperature-association model, temperature, relative laying date, age class and centred year enter location and scale. Female and year random intercepts enter location. There is no random year intercept in scale; the attempted additional term reached its lower optimisation bound with negligible likelihood improvement.
- **Egg-to-fledging success:** zero-inflated beta-binomial mixed model fitted by ML in glmmTMB. The same fixed predictors enter location and precision. Female and year random intercepts enter location; independent year random intercepts enter precision and extra-zero probability. The zero component has an intercept and year random effect.
- **Spring thermal conditions:** Gaussian location-scale ML model for 46 annual means (1980â€“2025), with centred year in mean and log residual SD.

The main reproductive models are additive. The site does not present exploratory age interaction analyses. Reproductive consistency is interpreted separately for each outcome distribution; the fledging precision measure is conditional on mean and clutch size and does not capture every source of full-mixture variability. Diagnostic limitations remain documented on the model pages.

## Files

| Location | Contents |
| --- | --- |
| `tables/final_models_20261019/` | Final COM-Poisson coefficients, predictions, centring constants and diagnostics |
| `tables/primary_ZIyear_ML_20261021/` | Final fledging-success estimates, prediction CIs, random-effect summaries and model checks; labelled tables also retain the simpler ML comparator |
| `tables/temporal_models_20261022/` | Temporal models without temperature: coefficients, prediction CIs, random SDs, diagnostics and simulation checks |
| `data_derived/` | Public environmental summaries and final model-ready breeding records in `analysis/` |
| `figures/` | Final published figure files |
| `R/` | Retained ML fitting/verification/export scripts from the original workflow |
| `scripts/` | Figure-generation code and fixed-design VIF check |
| `*.qmd`, `_quarto.yml`, `references.bib`, `styles.css` | Website source |
| `docs/` | Rendered website served by GitHub Pages |
| `publication_manifest.json` | SHA-256 checksums for numerical exports and figures |

## Rebuild the website

Install Quarto, then run from the repository root:

```sh
quarto render
```

The website uses the archived final exports. Rendering does not refit the reproductive models. GitHub Pages serves `main:/docs`.

## Scope of reproducibility and data access

The final curated model-ready records are public in `data_derived/analysis/`, with model-only female identifiers preserving repeated observations across both traits. See its README for the variable dictionary and sample-specific centring. The larger preparation CSV is not interchangeable with a final model sample.

```sh
Rscript R/95_validate_public_model_data.R
Rscript R/94_reproduce_from_public_data.R
```

The runner starts from CSV records and public numerical starting values, refits the main reproductive models, and writes regenerated objects and numerical exports separately from the published results. Historical fitted objects are not needed. R with glmmTMB >= 1.1.14, Matrix and MASS is required. See the website's Data preparation page for stages, output paths and numerical-validation limits. Original software versions remain recorded in the historical session-info files. Raw field-data cleaning is outside this model-ready release; annual phenology-window inputs are separately public.

## Publication integrity

The 22 October 2026 temporal revision replaces the displayed year curves with estimates from separately refitted models without temperature in either location or scale. It does not change the original temperature-association fits or the temperature, laying-date and age estimates and figures. Numerical exports and figure files are tracked by SHA-256 in `publication_manifest.json`. Superseded temperature-adjusted time-curve exports are retained in the authors' private archive and Git history rather than presented as current temporal results. All figure bands are ML confidence intervals, not prediction intervals for individual broods.

## Recovered sliding-window workflow (26 October 2026)

`R/04_climwin_laying_date.R` was recovered from the original working project and copied unchanged. `cache/climwin/recovered_search_data.rds` is a data-only extract of its aggregate search caches, including all 500 cached randomised searches. `R/04a_audit_recovered_climwin.R` checks these outputs and reconstructs complete-week exposures from public weather data without repeating a search. It reproduces the reported held-out slope and predictive RMSE. The historical `singlewin` output used a different exposure series; that output is retained unchanged under `tables/climwin_recovered/archive/`, separately from the checked results. Source checksums and recovery session information accompany the exports.

Run `Rscript R/04a_audit_recovered_climwin.R` and `Rscript R/04b_plot_recovered_climwin.R` from the repository root to rebuild the recovered-output tables and diagnostic figure with base R. The original full-search script requires the private raw bird input named in the script and its package environment. The website reads saved outputs and never repeats the expensive search during rendering.

`R/03_predictor_diagnostics.R` computes fixed-design VIF and year-temperature correlations from the public curated model-ready CSVs, verifying their annual sample summaries against the final exports. Both model-ready records and aggregate diagnostics are public.

## Original CF source archive

The unchanged original 1980–2019 field-record file is provided in [data_raw/birds/CF_adult_birds_1980-2019.csv](data_raw/birds/CF_adult_birds_1980-2019.csv). See its [provenance documentation](data_raw/birds/README.md). It precedes the later breeding dataset used in preparation and is not a substitute for the final model-ready inputs in `data_derived/analysis/`.
