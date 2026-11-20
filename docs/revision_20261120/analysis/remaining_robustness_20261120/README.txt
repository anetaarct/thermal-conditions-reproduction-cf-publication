Additional robustness checks completed on the corrected samples, 20 November 2026

Validated: repeated-female pooled clutch model; pooled and within-between
fledging models; full-sample clutch pooled and within-between models with an
added independent year random intercept in scale. Their selected BFGS fits have
convergence code 0, positive-definite Hessians and scaled gradients below 0.005.
The initial pooled-clutch nlminb warning was resolved by an independent restart
(independent_restart_diagnostics.csv). Scale-year nlminb singular warnings
occurred as the variance approached zero; BFGS refinements reproduced the
primary likelihood and coefficients, with SD about 7e-7. Do not interpret a
non-zero extra variance or use a regular interior-parameter variance LRT.

Not validated: repeated-female clutch within-between model. Multiple starting
values, BFGS, nlminb and L-BFGS-B, including finite-difference gradients, were
attempted. An independent audit of a reconstructed nearby candidate parameter
point found a maximum absolute analytic gradient about 9.2e7 and an objective
change about 0.012 at the same fixed parameters. Stalled/affected attempts were
terminated. No estimate, AIC or confidence interval from that model is used.
This check is numerically inconclusive for the primary within-female SD result.
Absence of a validated fit is not evidence for or against the temperature effect.

The repeated samples contain 6169 clutch records/2500 females and 4564 fledging
records/1907 females. Original centring, relative date and thermal decomposition
are preserved. Do not compare AIC between restricted and full samples.

Per-task directories contain the actual numerical exports and failed-attempt
logs. Some logs end in a trailing parse error after COMPLETE because the driver
file was revised while an earlier Rscript process was streaming it. All five
reported selected model objects and exports were already saved; their actual
diagnostics and independent likelihood agreement were checked separately.
Initial scaling attempts were stopped before producing inferential results.
Only validated exports are inferential results; failed logs are audit records.

Reproduction: after the main corrected fitting runner, run
Rscript R/revision_20261119/fit_remaining_robustness.R
This fits the five retained specifications, with outputs under
remaining_robustness_20261120/all/. It explicitly skips the unvalidated clutch
within-between specification by default. The separate bounded/numerical scripts
document attempts on that specification and expect the per-task directory
layout recorded here. RUN_UNVALIDATED_CLUTCH_WB=true enables the corresponding
unvalidated attempt in the main robustness script. It may stall and must not be
treated as an accepted fit without independent validation.

The updated assembled public wrapper has not been executed end-to-end in a
further fresh checkout. This limitation is unchanged. No raw-data cleaning
workflow or private ring mapping is added by this update.
