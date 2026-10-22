# Warmer springs and reproductive performance in collared flycatchers

Publication materials for *Larger clutches but lower egg-to-fledging success in warmer springs in a migratory passerine*.

**Results website:** https://anetaarct.github.io/thermal-conditions-reproduction-cf-publication/

This repository starts with a clean publication history. It contains the final maximum-likelihood analyses reported in the manuscript, their exported results and the companion Quarto website. Working drafts and superseded analyses are kept separately in a private archive.

## Reported models

- **Clutch size:** mean-parameterised COM-Poisson mixed model fitted by ML in glmmTMB. Temperature, relative laying date, age class and centred year enter location and scale. Female and year random intercepts enter location. There is no random year intercept in scale; the attempted additional term reached its lower optimisation bound with negligible likelihood improvement.
- **Egg-to-fledging success:** zero-inflated beta-binomial mixed model fitted by ML in glmmTMB. The same fixed predictors enter location and precision. Female and year random intercepts enter location; independent year random intercepts enter precision and extra-zero probability. The zero component has an intercept and year random effect.
- **Spring thermal conditions:** Gaussian location-scale ML model for 46 annual means (1980–2025), with centred year in mean and log residual SD.

The main reproductive models are additive. The site does not present exploratory age interaction analyses. Reproductive consistency is interpreted separately for each outcome distribution; the fledging precision measure is conditional on mean and clutch size and does not capture every source of full-mixture variability. Diagnostic limitations remain documented on the model pages.

## Files

| Location | Contents |
| --- | --- |
| `tables/final_models_20261019/` | Final COM-Poisson coefficients, predictions, centring constants and diagnostics |
| `tables/primary_ZIyear_ML_20261021/` | Final fledging-success estimates, prediction CIs, random-effect summaries and model checks; labelled tables also retain the simpler ML comparator |
| `data_derived/` | Public aggregate environmental and annual summary data |
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

Individual-level field records and fitted model objects are not included in this public repository. Curated brood-level data can be made available to reviewers on request. The retained R scripts record the original ML workflow and may refer to the authors' computational paths, saved inputs and fitted objects. They are not a standalone, one-command refit from the public files alone. Full model refitting requires the original curated individual-level data, centring constants, fitted starting objects where used, and the recorded software environment. The figure script uses public prediction exports and requires Python with pandas, NumPy and Matplotlib. It is supplied as code; this migration did not rerun it or change any published figure.

Daily weather observations originate from the Swedish Meteorological and Hydrological Institute (SMHI). Field-data enquiries should be directed to the corresponding author listed in the manuscript.

## Publication integrity

During the move to this repository, numerical CSV exports, figure files and the retained R scripts were copied without modification and checked by SHA-256. Only repository/site links, website navigation metadata, the README and explanatory text were updated. No models were refitted. All figure bands and model-estimate intervals retain their original meaning; ML confidence intervals must not be interpreted as prediction intervals for individual broods.
