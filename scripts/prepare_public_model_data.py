"""Audit local curated CSVs against final exports and prepare model-only inputs."""
from pathlib import Path
import csv, json, hashlib, math
from collections import defaultdict
from statistics import mean, median, stdev

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parent
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def num(x): return float(x) if x not in ('','NA') else None
def close(a,b,tol=1e-9): return abs(float(a)-float(b))<=tol
clutch=read(SOURCE/'data_derived/clutch_size_model_data.csv')
repro=read(SOURCE/'data_derived/reproductive_model_data.csv')
old_fledging=read(SOURCE/'data_derived/fledging_model_data_1981_2019.csv')
fledging=[r.copy() for r in repro if int(r['year'])<=2019 and num(r['fledged']) is not None
          and r['age_class'] in ('YOUNG','OLD') and r['female_id'] and num(r['rel_LD']) is not None]
byyear=defaultdict(list)
for r in fledging: byyear[r['year']].append(r)
tc=mean(float(r['temp_mean']) for r in fledging)
yc=mean(int(r['year']) for r in fledging)
for rows in byyear.values():
    centre=median(float(r['laying_date']) for r in rows)
    for r in rows:
        r['rel_LD']=float(r['laying_date'])-centre
        r['temp_mean_c']=float(r['temp_mean'])-tc
        r['year_c']=int(r['year'])-yc
        r['year_f']=r['year']
for r in clutch: r['year']=r['year_f']
checks=[]
def check(name,passed,detail=''):
    checks.append(dict(check=name,passed=bool(passed),detail=detail))
    if not passed: raise AssertionError(name+': '+str(detail))
fields=['female_id','year_f','clutch','rel_LD','temp_mean_c','year_c','age_class']
archived=read(SOURCE/'glmmTMB_comparison_20261016/clutch_counts/data.csv')
check('Clutch rows match archived final model matrix',len(clutch)==len(archived) and all(
 all(a[k]==b[k] if k in ('female_id','year_f','age_class') else close(a[k],b[k]) for k in fields)
 for a,b in zip(clutch,archived)),f'{len(clutch)} rows')
check('Fledging CSV matches reconstructed final sample',len(fledging)==len(old_fledging) and all(
 all(a[k]==b[k] if k in ('female_id','year_f','age_class') else close(a[k],b[k])
     for k in fields+['fledged']) for a,b in zip(fledging,old_fledging)),f'{len(fledging)} rows')
reference=read(ROOT/'tables/temporal_models_20261022/observed_annual.csv')
meta=read(ROOT/'tables/within_between_20261026/sample_metadata.csv')
for trait,rows in [('clutch',clutch),('fledging',fledging)]:
    groups=defaultdict(list)
    for r in rows: groups[r['year']].append(r)
    expected=[r for r in reference if r['trait']==trait]
    for e in expected:
        g=groups[e['year']]
        response=[float(r['clutch']) if trait=='clutch' else float(r['fledged'])/float(r['clutch']) for r in g]
        for label,a,b in [('n',len(g),e['n']),('observed_mean',mean(response),e['observed_mean']),
            ('observed_SD',stdev(response),e['observed_SD']),
            ('older_fraction',mean(r['age_class']=='OLD' for r in g),e['older_fraction']),
            ('mean_relative_date',mean(float(r['rel_LD']) for r in g),e['mean_relative_date']),
            ('mean_temperature',mean(float(r['temp_mean_c']) for r in g),e['mean_temperature'])]:
            check(f'{trait} {e["year"]} {label}',close(a,b))
        if trait=='fledging': check(f'{trait} {e["year"]} zero_fraction',close(mean(float(r['fledged'])==0 for r in g),e['zero_fraction']))
    m=next(r for r in meta if r['trait']==trait)
    check(trait+' final counts',len(rows)==int(m['records']) and len(groups)==int(m['seasons']) and
          len({r['female_id'] for r in rows})==int(m['females']))
    check(trait+' centering',close(mean(float(r['temp_mean'])-float(r['temp_mean_c']) for r in rows),m['temperature_center'],1e-8)
          and close(mean(int(r['year'])-float(r['year_c']) for r in rows),m['year_center'],1e-8))

# Model-only identifiers preserve repeated observations and ordering across traits.
ids=sorted({r['female_id'] for r in clutch+fledging})
idmap={v:f'F{i:05d}' for i,v in enumerate(ids,1)}
out=ROOT/'data_derived/analysis';out.mkdir(exist_ok=True)
columns=['female_id','year','year_f','age_class','laying_date','rel_LD','temp_mean','temp_mean_c','year_c','clutch']
# Clutch loading CSV lacks absolute laying date; join it by source row without retaining that source key.
source_rows={r['source_row']:r for r in repro}
for r in clutch: r['laying_date']=source_rows[r['source_row']]['laying_date']
hashes={}
for trait,rows in [('clutch',clutch),('fledging',fledging)]:
    names=columns+(['fledged'] if trait=='fledging' else [])
    path=out/f'{trait}_model_data.csv'
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=names);w.writeheader()
        for r in rows:
            z={k:r[k] for k in names};z['female_id']=idmap[r['female_id']];w.writerow(z)
    hashes[str(path.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(path.read_bytes()).hexdigest()
with (out/'checksums.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['filename','md5','sha256']);w.writeheader()
    for path,h in hashes.items():
        p=ROOT/path;w.writerow(dict(filename=p.name,md5=hashlib.md5(p.read_bytes()).hexdigest(),sha256=h))
report={'status':'passed','checks':len(checks),'results':checks,'sha256':hashes,
        'source_note':'Final-model rows and predictors verified against archived clutch matrix, reconstructed fledging sample, final annual exports and sample metadata. Fitted objects unavailable locally; this is not a bytewise comparison with their model frames.',
        'broader_reproductive_csv':'Preparation input, not a final analysis sample. Fledging requires sample-specific date and covariate recentering.'}
(ROOT/'tables/model_data_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':report['status'],'checks':len(checks),'sha256':hashes},indent=2))
