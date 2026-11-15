"""Refresh display model inputs after a licensed-data analysis run."""
from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parents[2]
for name in ['endpoints','icu','classifier']:
    for source in (ROOT/'src'/name/'aggregate').glob('*.csv'):
        destination=ROOT/'data/aggregate'/name/source.name
        shutil.copy2(source,destination)
for name in ['all_icu_gbmtm_selected_fitted_profiles.csv','all_icu_gbmtm_candidate_metrics_k1_k6.csv']:
    shutil.copy2(ROOT/'results/trajectory/tables'/name,ROOT/'data/aggregate/trajectory'/name)

shutil.copy2(ROOT/'src/classifier/aggregate/classifier_manifest.json',ROOT/'data/aggregate/classifier/classifier_manifest.json')
for source in (ROOT/'results/trajectory/tables').glob('all_icu_gbmtm_parameter_sensitivity*.csv'):
    shutil.copy2(source,ROOT/'data/aggregate/trajectory'/source.name)
for source in (ROOT/'results/trajectory/tables').glob('all_icu_gbmtm_k*_multistart_stability.csv'):
    shutil.copy2(source,ROOT/'data/aggregate/trajectory'/source.name)
