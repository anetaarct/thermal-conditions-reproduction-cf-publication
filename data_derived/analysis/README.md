# Final model-ready breeding records

These two CSV files contain the curated observations used by the final reproductive models, rather than the larger preparation dataset. Run commands from the repository root.

* `clutch_model_data.csv`: 12,627 attempts, 8,894 females, 44 seasons (1981–2025, excluding 2020).
* `fledging_model_data.csv`: 9,790 attempts, 7,132 females, 39 seasons (1981–2019).

`female_id` is a model-only identifier, consistent across both files. Original ring labels and source-row identifiers are omitted. Repeated observations and grouping are preserved. Each row is one included breeding attempt. The identifiers are ordered consistently with the original model grouping levels.

| Column | Definition |
| --- | --- |
| female_id | Shared female grouping identifier |
| year, year_f | Calendar year; use year_f as a factor for random effects |
| age_class | YOUNG = yearling; OLD = at least two years old |
| laying_date | First-egg date in days relative to 1 May (0 = 1 May) |
| rel_LD | First-egg date minus the relevant sample's annual median |
| temp_mean | Mean spring air temperature, °C, 26 March–7 May |
| temp_mean_c | Temperature minus the trait-specific centring constant |
| year_c | Year minus the trait-specific centring constant |
| clutch | Initial egg count; binomial denominator for fledging success |
| fledged | Fledgling count, including complete failures; fledging file only |

Clutch year is centred at 2002.61363636 and temperature at 5.1325594288°C. Fledging year is centred at 2004.19785495 and temperature at 5.0985959807°C. Do not recenter using an unweighted mean of annual temperatures. The `rel_LD` values use the final trait-specific samples; the earlier preparation dataset used a broader median for fledging.

The author-side audit (`tables/model_data_audit.json`) records 543 checks, source comparisons and SHA-256 checksums. The public validation script checks the released records against the final annual exports and within–between sample metadata. The historical fitted objects were not available for a bytewise model-frame comparison; numerical refitting provides an additional check.

```sh
Rscript R/95_validate_public_model_data.R
Rscript R/94_reproduce_from_public_data.R
```

The complete run requires R, glmmTMB >= 1.1.14, Matrix and MASS. The publication's original software versions are recorded in the relevant session-info files. Optimisation and Monte Carlo intervals can differ slightly across versions and platforms. The runner uses published coefficients as numerical initial values, not fitted-model objects. Models are refitted by likelihood optimisation with numerical checks. Outputs and newly generated model checkpoints are kept separately from the publication exports.

The raw field-data cleaning history is not reproduced from these files. They start at the final curated samples. Sliding-window selection uses the public annual laying-date summaries, including the seasons outside these trait-specific samples.
