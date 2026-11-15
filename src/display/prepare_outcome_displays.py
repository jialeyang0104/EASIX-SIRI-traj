#!/usr/bin/env python3
"""Bind the approved 28-day/ICU endpoints to the current article displays."""
from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
W = HERE.parent
A = W / 'endpoints/aggregate'
ICU = W / 'icu/aggregate'
S = HERE.parents[1] / 'data/aggregate/display'
EXPOSURES = ['nee_any_24h', 'nee_high', 'fluid_kg_high', 'documented_mv_24h', 'documented_crrt_24h']
COLORS = {1: '#4C78A8', 2: '#59A14F', 3: '#E15759'}


def publish(frame, name):
    target = S / f'{name}.csv'
    temporary = target.with_suffix('.csv.tmp')
    frame.to_csv(temporary, index=False)
    temporary.replace(target)


def build_mortality():
    static = pd.read_csv(A / 'all_icu_cox_static_effects.csv')
    dynamic = pd.read_csv(A / 'all_icu_cox_gbmtm_trajectory_effects.csv')
    static = static.loc[static.analysis.eq('static_day1_24h_landmark') & static.model.eq('static_joint')
                        & static.horizon.eq('28-day') & static.term.isin(['log2_easix_z', 'log2_siri_z'])].copy()
    dynamic = dynamic.loc[dynamic.analysis.eq('day7_landmark') & dynamic.horizon.eq('28-day')
                          & dynamic.term.isin(['trajectory_class2', 'trajectory_class3'])].copy()
    static['source_path'] = 'src/endpoints/aggregate/all_icu_cox_static_effects.csv'
    dynamic['source_path'] = 'src/endpoints/aggregate/all_icu_cox_gbmtm_trajectory_effects.csv'
    original = pd.concat([static, dynamic], ignore_index=True)
    original['effect_measure'] = 'HR'
    icu = pd.read_csv(ICU / 'fine_gray_estimates.csv')
    icu = icu.loc[icu.model.isin(['static_joint', 'gbmtm_k3'])].copy()
    icu['horizon'] = 'ICU'
    icu['hazard_ratio'] = icu.effect
    icu['events'] = icu.icu_deaths
    icu['source_path'] = 'src/icu/aggregate/fine_gray_estimates.csv'
    result = pd.concat([original, icu], ignore_index=True)
    labels = {'log2_easix_z': 'EASIX', 'log2_siri_z': 'SIRI',
              'trajectory_class2': 'Class 2 vs 1', 'trajectory_class3': 'Class 3 vs 1'}
    result['display_label'] = result.term.map(labels)
    assert len(result) == 16 and result.display_label.notna().all()
    publish(result, 'figure_5_mortality_forest')


def build_subgroups():
    effects = pd.read_csv(A / 'all_icu_cox_gbmtm_subgroup_stratified_effects.csv')
    tests = pd.read_csv(A / 'all_icu_cox_gbmtm_subgroup_interactions.csv')
    effects = effects.loc[effects.analysis.eq('day7_landmark') & effects.horizon.eq('28-day')].copy()
    tests = tests.loc[tests.analysis.eq('day7_landmark') & tests.horizon.eq('28-day')].copy()
    effects['effect_measure'] = 'HR'
    tests['test_method'] = 'Likelihood-ratio'
    ie = pd.read_csv(ICU / 'subgroup_estimates.csv')
    it = pd.read_csv(ICU / 'subgroup_interaction_tests.csv')
    ie['horizon'] = 'ICU'
    ie['hazard_ratio'] = ie.effect
    ie['events'] = ie.icu_deaths
    it['horizon'] = 'ICU'
    it['events'] = it.icu_deaths
    it['test_method'] = 'Robust joint Wald'
    publish(pd.concat([effects, ie], ignore_index=True), 'figure_7_subgroup_effects')
    publish(pd.concat([tests, it], ignore_index=True), 'figure_7_subgroup_interactions')


def build_treatments():
    estimates = pd.read_csv(A / 'supp_treatment_forest_estimates.csv')
    tests = pd.read_csv(A / 'supp_treatment_interaction_tests.csv')
    support = pd.read_csv(A / 'supp_treatment_support_cells.csv')
    estimates = estimates.loc[estimates.horizon.eq(28) & estimates.exposure.isin(EXPOSURES)].copy()
    tests = tests.loc[tests.horizon.eq(28) & tests.exposure.isin(EXPOSURES)].copy()
    support = support.loc[support.horizon.eq(28) & support.exposure.isin(EXPOSURES)].copy()
    for frame in [estimates, tests, support]:
        frame['horizon'] = '28-day'
    estimates['effect_measure'] = 'HR'
    tests['effect_measure'] = 'HR'
    ie = pd.read_csv(ICU / 'treatment_estimates.csv')
    it = pd.read_csv(ICU / 'treatment_interaction_tests.csv')
    ic = pd.read_csv(ICU / 'treatment_support_cells.csv')
    ie = ie.loc[ie.exposure.isin(EXPOSURES)].copy()
    it = it.loc[it.exposure.isin(EXPOSURES)].copy()
    ic = ic.loc[ic.exposure.isin(EXPOSURES)].copy()
    ie['horizon'] = 'ICU'
    ie['hazard_ratio'] = ie.effect
    ie['model_status'] = ie.status
    ie['model_n'] = ie.n
    ie['model_events'] = ie.icu_deaths
    ie['global_interaction_p'] = ie.interaction_p
    ie['source_path'] = 'src/icu/aggregate/treatment_estimates.csv'
    ie['encoding'] = 'cutoff_binary'
    it['horizon'] = 'ICU'
    it['model_n'] = it.n
    it['model_events'] = it.icu_deaths
    it['global_interaction_p'] = it.interaction_p
    it['interaction_df'] = it.df
    it['interaction_chisq'] = it.wald_chisq
    it['effect_measure'] = 'sHR'
    it['source_path'] = 'src/icu/aggregate/treatment_interaction_tests.csv'
    it['encoding'] = 'cutoff_binary'
    ic['horizon'] = 'ICU'
    for original, new, name in [(estimates, ie, 'supp_treatment_forest_estimates'),
                                 (tests, it, 'supp_treatment_interaction_tests'),
                                 (support, ic, 'supp_treatment_support_cells')]:
        result = pd.concat([original, new], ignore_index=True)
        result['exposure_order'] = result.exposure.map({v: i for i, v in enumerate(EXPOSURES)})
        if 'trajectory_class' in result:
            result['class_color'] = result.trajectory_class.astype(int).map(COLORS)
        assert len(result) == {'supp_treatment_forest_estimates': 60,
                               'supp_treatment_interaction_tests': 20,
                               'supp_treatment_support_cells': 120}[name]
        publish(result, name)
    publish(pd.read_csv(S / 'supp_treatment_support_cells.csv'), 'supp_treatment_cutpoint_support')


def build_flow_counts():
    target = S / 'figure_1_cohort_flow.json'
    config = json.loads(target.read_text())
    manifest = pd.read_csv(ICU / 'model_manifest.csv')
    models = manifest.loc[manifest.role.eq('primary_association') & manifest.estimator.eq('Fine-Gray')]
    records = []
    for cohort, db in zip(config['cohorts'], ['mimiciv', 'mimiciii']):
        m = models.loc[models.dataset.eq(db)]
        cohort['analysis']['static_icu_models'] = {
            name: int(m.loc[m.model.eq('static_' + name), 'n'].iloc[0])
            for name in ['easix', 'siri', 'joint']}
        cohort['analysis']['trajectory_icu_day7'] = int(m.loc[m.model.eq('gbmtm_k3'), 'n'].iloc[0])
        for name, n in cohort['analysis']['static_icu_models'].items():
            records.append(dict(dataset=db, analysis='static_day1_24h_landmark', model=name,
                                outcome='ICU', n=n))
        records.append(dict(dataset=db, analysis='day7_landmark', model='gbmtm_k3', outcome='ICU',
                            n=cohort['analysis']['trajectory_icu_day7']))
    temporary = target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(config, ensure_ascii=False, indent=2) + '\n')
    temporary.replace(target)
    publish(pd.DataFrame(records), 'figure_1_icu_analysis_populations')


def main():
    build_mortality()
    build_subgroups()
    build_treatments()
    build_flow_counts()
    print('Bound 28-day Cox and ICU Fine–Gray displays: 16 main effects, 60 treatment effects, 20 treatment tests.')


if __name__ == '__main__':
    main()
