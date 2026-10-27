"""Build a code-and-output chapter after the within-between models have passed checks."""
from pathlib import Path
import csv

root = Path(__file__).resolve().parent.parent
folder = root / "tables/within_between_20261026"

def rows(name):
    with (folder / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def table(data, columns):
    def fmt(value):
        try:
            return f"{float(value):.5g}"
        except (TypeError, ValueError):
            return str(value)
    labels = {"females_varied_temperature": "Females with temperature variation",
              "records_in_varied_temperature_females": "Their records"}
    return ("::: {.table-responsive}\n\n| "
            + " | ".join(labels.get(c, c.replace("_", " ")) for c in columns)
            + " |\n| " + " | ".join("---" for _ in columns) + " |\n"
            + "".join("| " + " | ".join(fmt(row[c]) for c in columns) + " |\n" for row in data)
            + "\n:::\n\n")

def panel(label, body, opened=False):
    return (f'<details class="analysis-code"{" open" if opened else ""}>\n'
            f"<summary>{label}</summary>\n\n{body.strip()}\n\n</details>\n\n")

def include(script, title):
    source = root / script
    target = root / "_includes" / f"{source.stem}.qmd"
    text = source.read_text(encoding="utf-8-sig").rstrip()
    target.write_text(f'```{{.r filename="{script}"}}\n{text}\n```\n', encoding="utf-8")
    return panel(title, "{{< include _includes/" + target.name + " >}}")

fixed = rows("fixed_effects.csv")
metadata = rows("sample_metadata.csv")
comparison = rows("model_comparison.csv")
contrasts = rows("within_between_contrasts.csv")
checks = rows("conditional_simulation_checks.csv")
diagnostics = rows("optimizer_comparison.csv")
assert {r["trait"] for r in metadata} == {"clutch", "fledging"}
assert all(r["convergence"] == "0" and r["pdHess"] == "TRUE"
           and float(r["scaled_gradient"]) < .005 for r in diagnostics)
temperature = [r for r in fixed if r["model"] == "within_between"
               and r["term"] in ("temp_within", "temp_between")]
assert len(temperature) == 8

text = '''---
title: "Temperature within and between females"
number-sections: false
execute:
  enabled: false
---

## Data preparation

This second temperature formulation separates temperature associations within the same female from differences between females. Temperature is split in both location and scale; continuous year, laying date, age and the existing random-effect structures are retained.

'''
text += panel("Show temperature decomposition", '''```r
# d is the trait-specific analysis dataset; female means use its observed records.
centre <- unique(d$temp_mean - d$temp_mean_c)
d$temp_between <- ave(d$temp_mean, d$female_id, FUN = mean) - centre
d$temp_within <- d$temp_mean_c - d$temp_between
stopifnot(max(abs(d$temp_within + d$temp_between - d$temp_mean_c)) < 1e-8)
```''', True)
text += table(metadata, ["trait", "records", "females", "seasons", "females_varied_temperature",
                         "records_in_varied_temperature_females"])
text += ("Single-record females have zero within deviation. The between predictor describes the seasons in which a female was observed, not a stable thermal preference.\n\n"
         "## Fit the within-between models\n\n"
         "Clutch size uses COM-Poisson; fledging success uses zero-inflated beta-binomial. The extra-zero equation retains an intercept and year random effect, without temperature slopes.\n\n")
text += include("R/90_fit_temperature_within_between.R", "Show complete fitting and optimiser-check code")
text += "## Read saved results\n\n"
text += panel("Show cached-output loading", '''```r
results <- readRDS("tables/within_between_20261026/aggregate_results.rds")
results$fixed_effects
results$contrasts
```''', True)
text += table(temperature, ["trait", "component", "term", "estimate", "lower", "upper"])
text += ("Location (`cond`): log clutch mean or logit beta-binomial mean. Scale (`disp`): log COM-Poisson inverse shape or log beta-binomial precision. Response-scale consistency curves below account for both clutch components.\n\n"
         "[All fixed effects](tables/within_between_20261026/fixed_effects.csv) · "
         "[Aggregate RDS output](tables/within_between_20261026/aggregate_results.rds)\n\n"
         "## Within versus between: direct tests\n\n")
text += table(contrasts, ["trait", "component", "estimate", "lower", "upper", "p"])
text += ("Each contrast is **within minus between** on its component's link scale. A significant coefficient alongside a non-significant coefficient is not itself evidence that the effects differ.\n\n"
         "## Comparison with the single-temperature model\n\n")
text += table(comparison, ["trait", "baseline_AIC", "within_between_AIC", "delta_AIC", "LR", "df", "p"])
text += ("Both members of each comparison were fitted under R 4.6.1 / glmmTMB 1.1.15.2. The baseline likelihoods were checked against the original final fits. Splitting temperature adds two fixed coefficients in total; the nested comparison tests equality of within and between effects in both components.\n\n"
         "## Prediction export and checks\n\n")
text += include("R/91_export_temperature_within_between.R", "Show prediction export and conditional-simulation code")
text += panel("Numerical convergence and optimiser checks",
              table(diagnostics, ["trait", "model", "optimizer", "convergence", "pdHess", "scaled_gradient", "logLik"]))
text += panel("Conditional simulation checks", table(checks, ["trait", "statistic", "observed", "predicted", "lower", "upper"]))
flagged = [f'{r["trait"]}: `{r["statistic"]}`' for r in checks
           if not float(r["lower"]) <= float(r["observed"]) <= float(r["upper"])]
if flagged:
    text += ("Observed statistics outside the conditional 95% simulation intervals: "
             + "; ".join(flagged) + ". These discrepancies remain visible in the diagnostics.\n\n")
text += panel("Fixed-design VIF", table(rows("fixed_design_vif.csv"), ["trait", "predictor", "VIF"]))
text += ("Simulation intervals use 300 plug-in datasets at fitted parameters and group modes; they are not cross-validation or parameter-uncertainty intervals. Full model and coefficient outputs are retained even when evidence is uncertain.\n\n"
         "## Figures\n\n"
         "Curves hold year and relative laying date at zero, age at yearling and all random effects at zero. Bands are pointwise 95% CIs from 4,000 joint coefficient draws. In the within panels, the female's mean temperature is held at the sample centre; in between panels, within deviation is zero.\n\n"
         "![Additional clutch-size model: within and between mean responses and conditional SD.](figures/Fig_within_between_clutch.png)\n\n"
         "[Clutch PDF](figures/Fig_within_between_clutch.pdf)\n\n"
         "![Additional fledging-success model: within and between mean responses and component precision.](figures/Fig_within_between_fledging.png)\n\n"
         "[Fledging PDF](figures/Fig_within_between_fledging.pdf) · "
         "[Figure code — final chapter](figure_code.qmd#temperature-within-between)\n\n"
         "## Interpretation\n\n"
         "These are adjusted associations in a within-between model. Lower clutch SD and higher fledging precision indicate greater consistency; fledging precision is conditional on component mean and clutch size. Between-female temperature differences may reflect sampling years and observation histories, so the decomposition does not establish causation or individual thermal preferences.\n\n"
         "[Session information](tables/within_between_20261026/session_info.txt) · "
         "[Export session](tables/within_between_20261026/export_session_info.txt) · "
         "[Reference effects per degree](tables/within_between_20261026/reference_effects_per_degree.csv)\n")
conclusions = []
for trait, label in [("clutch", "Clutch size"), ("fledging", "Fledging success")]:
    for component, name in [("cond", "location"), ("disp", "scale")]:
        r = next(r for r in contrasts if r["trait"] == trait and r["component"] == component)
        p = float(r["p"])
        wording = "evidence of a within/between difference" if p < .05 else "no clear within/between difference"
        conclusions.append(f'- **{label}, {name}:** {wording} (direct contrast p = {p:.3g}).')
effects = rows("reference_effects_per_degree.csv")
for r in effects:
    r["Mean change (95% CI)"] = (f'{float(r["mean_change_per_C"]):.3g} '
        f'({float(r["mean_lower"]):.3g}, {float(r["mean_upper"]):.3g})')
    r["Consistency ratio (95% CI)"] = (f'{float(r["consistency_ratio_per_C"]):.3f} '
        f'({float(r["ratio_lower"]):.3f}, {float(r["ratio_upper"]):.3f})')
effect_table = table(effects, ["trait", "predictor", "Mean change (95% CI)",
                              "consistency_measure", "Consistency ratio (95% CI)"])
text = text.replace("## Interpretation\n\n", "## Main results\n\n" + "\n".join(conclusions)
                    + "\n\nAt the reference conditions above, for a +1°C change: mean differences are in eggs for clutch size and proportions for fledging success. Consistency ratios compare conditional SD or component precision after versus before that change.\n\n"
                    + effect_table + "## Interpretation\n\n")
(root / "within_between.qmd").write_text(text, encoding="utf-8")
include("R/92_plot_temperature_within_between.R", "")
print("Built additional within/between model chapter from checked outputs.")
