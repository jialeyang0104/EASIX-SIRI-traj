"""Public tree contracts; no licensed data or model fitting required."""
from pathlib import Path
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
print('PASS Python syntax, portable paths, aggregate column audit, synthetic ICU endpoints, cohort conservation')
