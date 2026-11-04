"""Build visible, source-linked review summaries from public aggregate CSVs."""
from pathlib import Path
import csv
import re
import json

ROOT = Path(__file__).resolve().parents[1]
WB = 'tables/within_between_20261026/'
TC = 'tables/temporal_models_20261022/'
CW = 'tables/climwin_recovered/'

def table(path, columns=None, predicate=None, limit=None):
    with (ROOT / path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        columns = columns or reader.fieldnames
        rows = [r for r in reader if predicate is None or predicate(r)]
    if limit:
        rows = rows[:limit]
    def fmt(v):
        try:
            return format(float(v), '.6g') if v else '—'
        except ValueError:
            return v.replace('|', '\\|')
    lines = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---'] * len(columns)) + ' |']
    lines += ['| ' + ' | '.join(fmt(r[c]) for c in columns) + ' |' for r in rows]
    return '\n'.join(lines) + f'\n\n[Download complete results (CSV)]({path})'

def summary(kind):
    if kind == 'sample':
        return ('**Sample check.** Clutch size: 12,627 attempts, 8,894 females, 44 seasons. '
                'Egg-to-fledging success: 9,790 attempts, 7,132 females, 39 seasons. '
                'Public model-ready CSVs are available in data_derived/analysis; these counts are the saved analysis-sample metadata.\n\n'
                + table(WB+'sample_metadata.csv', ['trait','records','females','seasons','females_varied_temperature']))
    if kind == 'weather':
        return ('**Input check.** Public inputs contain 16,802 daily records, 46 annual temperatures and 45 annual laying-date summaries. '
                'The operational exposure is 26 March–7 May, 43 days each year. '
                'Run the code to check your downloaded files; this summary is not a new execution log.')
    if kind == 'window':
        return ('**Saved window-search results.** The selected complete-week search window is 6–0 weeks (49 days), '
                'distinct from the 43-day reproductive exposure. These are the archived results; new randomisations may differ.\n\n'
                + table(CW+'selected_windows_audit.csv'))
    if kind == 'validation':
        return ('**Saved validation.** The held-out association remains negative; prediction uses the training coefficients. '
                'This validates the phenology window, not the reproductive models.\n\n'
                + table(CW+'temporal_validation_audit.csv'))
    if kind == 'random':
        return ('**Saved randomisation check.** None of 500 cached searches exceeded the observed improvement; '
                'the plus-one empirical probability is 1/501 = 0.001996.\n\n'
                '[Download all 500 saved searches](tables/climwin_recovered/randomised_windows.csv).')
    if kind == 'predictors':
        return ('**Predictor checks.** Saved VIFs are calculated from the curated brood records; '
                'annual correlations alone cannot recreate brood-level VIFs.\n\n'
                + table('tables/predictor_diagnostics/fixed_design_vif.csv'))
    if kind == 'environment':
        return ('**Environmental result.** Mean spring temperature increased by 0.0439°C/year '
                '(95% CI 0.0248–0.0630). There was no clear change in residual interannual variability. '
                'The browser code reports time effects per decade; multiply annual slopes by ten when comparing. '
                'This is an environmental model, not a reproductive model.')
    if kind in ('clutch', 'fledging'):
        interpretation = ('Warmer springs were associated with larger clutches. The single-temperature consistency association remains uncertain. '
                          'COM-Poisson dispersion is log inverse shape, not log residual SD.' if kind == 'clutch' else
                          'Warmer springs were associated with lower egg-to-fledging success. Temperature associations with precision remain uncertain. '
                          'Precision describes the beta-binomial component; it is not full-mixture variance.')
        return ('**Saved model summary.** '+interpretation+' Coefficients below are on the link scale; '
                'these are exported estimates, not a fresh model fit. `cond` = location, `disp` = dispersion, `zi` = extra-zero probability.\n\n'
                + table(WB+'fixed_effects.csv', ['component','term','estimate','SE','lower','upper'],
                        lambda r: r['trait']==kind and r['model']=='baseline'))
    if kind in ('clutch_export','fledging_export'):
        if kind == 'clutch_export':
            path = 'tables/final_models_20261019/COM_temperature_mean.csv'
        else:
            path = 'tables/primary_ZIyear_ML_20261021/ZIyear_temperature_success.csv'
        return ('**Saved prediction output (first five rows).** `estimate` is the reference-condition expected outcome; '
                '`CI_low` and `CI_high` are pointwise 95% confidence limits. These are predictions, not model coefficients. '
                'Other predictors and random effects are held at the reference values specified in the chapter.\n\n'
                + table(path, limit=5))
    if kind == 'temporal':
        return ('**Saved calendar-year results, temperature omitted.** These are reference-condition changes from separate '
                'models without temperature, not the temperature-adjusted year coefficients. Mean changes are eggs for clutch size '
                'and proportions for fledging success.\n\n'
                + table(TC+'temporal_changes.csv', ['trait','first_year','last_year','mean_change','mean_change_lower','mean_change_upper']))
    if kind == 'within':
        return ('**Saved within–between comparison.** Splitting temperature improved clutch-model fit; '
                'it did not improve fledging-model fit. The direct clutch dispersion contrast supports a difference between levels; '
                'between-female exposure reflects observation history.\n\n'
                + table(WB+'model_comparison.csv')+'\n\n'
                + table(WB+'within_between_contrasts.csv'))
    if kind == 'joint':
        return ('**Temperature-adjusted calendar-year results.** Both formulations retain temperature. '
                'Mean changes are eggs or success proportions; consistency ratios describe clutch SD or fledging component precision.\n\n'
                + table('tables/joint_calendar_year_20261027/contrasts.csv',
                        ['trait','model','mean_change','mean_lower','mean_upper','consistency_ratio','ratio_lower','ratio_upper']))
    if kind == 'diagnostics':
        return ('**Saved numerical checks.** Convergence and a positive-definite Hessian assess numerical fit; '
                'they do not establish biological validity or out-of-sample prediction. Conditional simulations still underestimate '
                'overall reproductive variance.\n\n'+table(WB+'optimizer_comparison.csv',
                    ['trait','model','optimizer','AIC','convergence','pdHess','scaled_gradient']))
    return ('**Figure output.** This code draws saved estimates and pointwise 95% confidence intervals; '
            'it does not refit the models. The figure and its downloadable PDF are shown below. '
            'See the corresponding analysis chapter for coefficients and diagnostic limitations.')

def block(kind):
    return '\n\n::: {.callout-note .review-summary appearance="simple"}\n## Saved results for this step\n\n' + summary(kind) + '\n:::\n\n'

def choose(page, content, position, full):
    if 'lab-' in content:
        number = re.search(r'lab-(\d+)',content).group(1)
        return {'01':'weather','02':'window','03':'window','04':'environment','05':'predictors',
                '06':'clutch','07':'fledging','08':'within','09':'joint','10':'joint',
                '11':'figure','12':'figure','13':'temporal'}[number]
    if page == 'figure_code': return 'figure'
    if page == 'workflow': return 'weather'
    if page == 'data_preparation': return 'sample'
    if page in ('clutch_size','fledging_success'):
        trait = 'clutch' if page == 'clutch_size' else 'fledging'
        if 'Show R code' in content and 'readRDS' in content: return 'sample'
        if 'prediction-export' in content: return trait+'_export'
    if page == 'climate_window':
        if '04a_' in content or 'temporal_validation' in content: return 'validation'
        if 'p_plus_one' in content: return 'random'
        if 'nrow(annual)' in content: return 'weather'
        return 'window'
    if page == 'within_between':
        if '93_export' in content: return 'joint'
        if 'diagnostic' in content: return 'diagnostics'
        return 'within'
    if page == 'temporal_trends':
        return 'sample' if 'Show data-loading' in content else 'temporal'
    return page if page in ('clutch','fledging') else {'clutch_size':'clutch','fledging_success':'fledging'}.get(page,'environment')

counts = {}
for page in ('workflow','climate_window','data_preparation','model_structure','clutch_size','fledging_success','temporal_trends','within_between','figure_code'):
    path = ROOT / (page+'.qmd')
    original = path.read_text(encoding='utf-8')
    original = re.sub(r'\n\n::: \{\.callout-note \.review-summary[^\n]*\}\n.*?\n:::\n\n', '\n\n', original, flags=re.S)
    count = [0]
    # Match a complete code details block, or a browser-code include outside one.
    pattern = r'<details\b[^>]*>.*?</details>|\{\{< include _includes/lab-[^\n]+? >\}\}|```(?:r|bash)\n.*?\n```'
    def replace(m):
        text = m.group(0)
        if text.startswith('```') and ('install.packages' in text or 'R/94_reproduce_from_public_data.R' in text or 'R/95_validate_public_model_data.R' in text):
            return text  # Setup/entry-point commands, not an additional analysis.
        if text.startswith('<details') and ('```' not in text and '{{< include' not in text):
            # Keep saved output visible instead of hiding it beneath another click.
            return text.replace('class="analysis-code"','class="analysis-code" open',1) if 'open' not in text.split('>')[0] else text
        count[0] += 1
        return text + block(choose(page,text,m.start(),original))
    updated = re.sub(pattern,replace,original,flags=re.S)
    # Separate headings from preceding download paragraphs for correct Markdown rendering.
    updated = re.sub(r'([^\n])\n(## )',r'\1\n\n\2',updated)
    path.write_text(updated,encoding='utf-8')
    counts[page] = count[0]

(ROOT/'tables/review_summary_manifest.json').write_text(json.dumps(counts,indent=2)+'\n')
print(json.dumps(counts,indent=2))
