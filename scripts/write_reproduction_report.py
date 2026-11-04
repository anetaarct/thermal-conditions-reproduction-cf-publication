"""Summarise completed numerical checks without republishing fitted objects."""
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parents[1]
def read(p):
    with (ROOT/p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def comparison(a,b,keys):
    old=read(a);new=read(b);key=lambda r:tuple(r[k] for k in keys)
    new={key(r):r for r in new}
    old=[r for r in old if r.get('model','without_temperature')=='without_temperature'] if 'temporal' in a else old
    values=[abs(float(r['estimate'])-float(new[key(r)]['estimate'])) for r in old]
    assert len(values)==len(old) and max(values)<.001
    return dict(coefficients=len(values),max_absolute_difference=max(values),tolerance=.001,passed=True)
audit=json.loads((ROOT/'tables/model_data_audit.json').read_text())
report={'status':'passed','input_checks':audit['checks'],'data_sha256':audit['sha256'],
        'baseline_and_within_between':comparison('tables/within_between_20261026/fixed_effects.csv',
            'tables/within_between_20261026_reproduced/fixed_effects.csv',['trait','model','component','term']),
        'temporal_without_temperature':comparison('tables/temporal_models_20261022/fixed_effects.csv',
            'models/reproduced_temporal/fixed_effects.csv',['trait','model','component','term']),
        'window_search':read('tables/reproduced/selected_window.csv'),
        'environment_model':read('tables/reproduced/environment_coefficients.csv'),
        'poisson_comparator':read('tables/reproduced/poisson_diagnostics.csv'),
        'simulation_scope':'Conditional plug-in checks; regenerated Monte Carlo intervals may vary across numerical environments. Not held-out validation.',
        'scope':'Analyses start at final curated samples. Raw field cleaning is outside this release. No historical fitted-model object is required.'}
assert len(read('tables/reproduced/fresh_randomised_searches.csv'))==500
report['fresh_randomised_searches']=500
(ROOT/'tables/reproduction_validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:report[k] for k in ('status','baseline_and_within_between','temporal_without_temperature')},indent=2))
