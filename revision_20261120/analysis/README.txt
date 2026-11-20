Supplementary data for the curation and sensitivity revision of 19 November 2026

This package contains corrected anonymous model-ready samples, numerical model
exports, exact-age subsets and scripts. It accompanies the revised manuscript
and is distinct from the earlier GitHub release.

The corrected samples exclude the entire recorded history of any female with
an experimental flag anywhere in CF21980-2025.csv. All records from ambiguous
duplicate female-season groups are also excluded. Relative laying dates and
centred covariates are recomputed within each corrected trait-specific sample.
See data_derived/analysis_clean_20261119/curation_audit.csv for counts.

Reproduction from model-ready data

1. Download a fresh copy of the original project repository at commit
   63b46bbf74e527d195729fadfeddd741441b2ff1.
2. Extract this package's contents into that repository root. Use a fresh copy
   with no previously generated model checkpoints.
3. Install R, glmmTMB >= 1.1.14, MASS and Matrix.
4. From the repository root, run
   Rscript run_revision_from_public_data.R

The runner uses the archived original coefficients as numerical starting values,
not as numerical targets for changed samples. Corrected annual sample checks
replace the old-sample checks. It generates fresh fitted objects, both temperature
formulations, temporal models without temperature, prediction exports, six age
sensitivity fits and independent optimiser checks, extra-zero comparisons,
fixed-window sensitivity fits and revised figure PDFs.
Fitted-model objects and original ring mappings are not distributed here.

Scope: this runner reproduces model fitting and figure exports from the supplied
model-ready records. It does not reconstruct the raw-field-data cleaning or the
exact-age evidence audit from the source archive. The exclusion-count audit and
the age-reconstruction assumptions are provided, but a complete public raw-data-
to-model-ready cleaning workflow is outside this package, as in the original
model-ready release. Do not describe the wrapper as a fully tested raw-to-results
pipeline. Its fitting steps were run in the separate Totoro revision checkout;
the assembled public wrapper was syntax-checked, and its figure scripts were
executed against the supplied exports, but it has not been rerun end-to-end as
an extracted archive in a further fresh checkout.

To redraw the eight revised reproductive figures from the supplied numerical
exports without refitting, run these from the repository root after extraction:
Rscript R/revision_20261119/plot_revised_results.R
Rscript R/revision_20261119/plot_window_sensitivity.R

Exact-age records

age_sensitivity_20261119/*_known_age_model_data.csv contains the age-known
subsets. Exact age was reconstructed from recorded hatch years, yearling records
and explicit integer AGEF codes from 1 to 10. Half-integer, plus and uncertain
codes were not interpreted as exact ages. Conflicting and impossible histories
were excluded. Integer-code interpretation is an explicit audit assumption.
Age classes, linear exact age and quadratic exact age were compared on identical
subsets, preserving the full corrected sample's thermal decomposition. These
are sensitivity analyses of conditional temperature associations, not estimates
of individual ageing trajectories. The central clutch-SD association is
uncertain in the known-age subset even with categorical age.

The manuscript reports validated age and extra-zero tables. Files documenting
initial failed fits are retained for audit and are not inferential results.
Reference curves use 4,000 joint coefficient draws. Temporal, explicit SD, age
and failure contrasts use 6,000. Plug-in simulations omit parameter uncertainty
and do not establish held-out prediction.

The original window-search annual summaries and environmental temperature
analysis belong to the archived project. New fixed-window sensitivity models
compare 26 March-7 May, 1 April-7 May, 1-30 April and 15 April-7 May on identical
corrected samples. The positive mean clutch-size association persists; success
and clutch-SD associations are less uniformly supported across exposures.
All retained fits have numerical diagnostics and independent optimiser checks.
Failed checks are audit records, not inferential fits. Exact-age and temperature-
free temporal results are supplementary analyses. The latter provide temporal
context without identifying the causal contribution of warming.


Additional robustness update — 20 November 2026

Five additional specifications were numerically validated; the repeated-female
clutch within-between model was attempted but not validated. See
remaining_robustness_20261120/README.txt and SM S1.14, Tables S19–S21.
The added scale-year variance approached zero. No robustness claim is made
for the failed repeated-female clutch within-between specification.
