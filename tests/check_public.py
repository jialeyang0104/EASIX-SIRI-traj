"""Public tree contracts; no licensed data or model fitting required."""
from pathlib import Path
from collections import defaultdict
import ast,csv,importlib.util,sys
ROOT=Path(__file__).resolve().parents[1]
for path in (ROOT/'src').rglob('*.py'):
    ast.parse(path.read_text(),filename=str(path))
    assert '/Users/' not in path.read_text(),path
for path in (ROOT/'data/aggregate').rglob('*.csv'):
    with path.open() as handle:columns=next(csv.reader(handle))
    assert not set(columns)&{'subject_id','hadm_id','icu_stay_id','icustay_id','stay_id','intime','outtime','admittime','deathtime','dod'},path
sys.path.insert(0,str(ROOT/'src/display'))
spec=importlib.util.spec_from_file_location('icu',ROOT/'src/icu/prepare_analysis.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
module.endpoint_tests()
from cohort_support import validate_config
import json
validate_config(json.loads((ROOT/'data/aggregate/display/figure_1_cohort_flow.json').read_text()))
display = ROOT / 'data/aggregate/display'
def rows(stem):
    with (display / f'{stem}.csv').open() as f:
        return list(csv.DictReader(f))

def key(row):
    return row['dataset'], row['outcome'], row['exposure']

def valid_ratio(row):
    return 0 < float(row['ci_low']) <= float(row['ratio']) <= float(row['ci_high']) < float('inf')

cohorts = {'mimiciv', 'mimiciii'}
outcomes = {'28-day', 'ICU'}
fixed_exposures = {
    'nee_any_24h': None, 'nee_dose': .30, 'net_fluid': 35, 'documented_crrt_24h': None
}
fixed = rows('figure_s6_fixed')
assert len(fixed) == 48
assert {(*key(r),r['trajectory_class']) for r in fixed} == {
    (db, end, t, cls) for db in cohorts for end in outcomes for t in fixed_exposures for cls in ('1','2','3')
}
for r in fixed:
    assert (float(r['cutoff']) if r['cutoff'] else None) == fixed_exposures[r['exposure']]
    if (r['dataset'],r['exposure'],r['trajectory_class']) == ('mimiciii','documented_crrt_24h','1'):
        assert r['status'] == 'not_estimable' and r['high_n'] == '0' and not r['ratio']
    else:
        assert r['status'] == 'estimated' and valid_ratio(r)
    assert 0 <= float(r['interaction_p']) <= 1

tests = rows('figure_s7_four_treatment_tests')
tests_by_cell = {key(r):r for r in tests}
prob_exposures = {'exposure_any_nee': None, 'exposure_nee_dose': .30, 'exposure_net_fluid':35, 'exposure_crrt':None}
assert len(tests) == 16 and set(tests_by_cell) == {
    (db, end, t) for db in cohorts for end in outcomes for t in prob_exposures
}
p_fields = ('interaction_p','slope_difference_p')
for r in tests:
    assert r['status'] == 'estimated' and r['interaction_df'] == '2'
    assert (float(r['cutoff']) if r['cutoff'] else None) == prob_exposures[r['exposure']]
    assert all(0 <= float(r[p]) <= 1 for p in p_fields)
    assert int(r['low_n']) + int(r['high_n']) == int(r['n'])
    assert int(r['low_events']) + int(r['high_events']) == int(r['events'])

grids = defaultdict(list)
for r in rows('figure_s7_four_treatment_curves'):
    t = tests_by_cell[key(r)]
    assert all(r[p] == t[p] for p in p_fields)
    assert valid_ratio(r) and 0 <= float(r['probability']) <= 1
    grids[(*key(r),r['predicted_class'])].append(float(r['probability']))
assert set(grids) == {(*k,c) for k in tests_by_cell for c in ('2','3')}
assert all(len(g) == 101 and g == sorted(set(g)) for g in grids.values())

histograms = defaultdict(list)
for r in rows('figure_s7_four_treatment_histograms'):
    histograms[(*key(r),r['predicted_class'],r['exposure_arm'])].append(r)
assert set(histograms) == {(*g,a) for g in grids for a in ('0','1')}
for group, h in histograms.items():
    n = int(tests_by_cell[group[:3]]['low_n' if group[-1] == '0' else 'high_n'])
    assert len(h) == 20 and sum(int(r['count']) for r in h) == n
    for i,r in enumerate(sorted(h,key=lambda r:float(r['bin_lower']))):
        assert abs(float(r['bin_lower'])-i/20) < 1e-12
        assert abs(float(r['bin_upper'])-(i+1)/20) < 1e-12

print('PASS Python syntax, portable paths, aggregate column audit, synthetic ICU endpoints, cohort conservation and S6/S7 display sources')
