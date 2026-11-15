#!/usr/bin/env python3
"""Prepare the independent ICU-death endpoint from current local source caches."""
from pathlib import Path
from datetime import datetime, timezone
import importlib.util
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"display"))
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
WORK = HERE.parent
SOURCE = WORK / 'endpoints'
PRIVATE = HERE / 'data_private'
AGGREGATE = HERE / 'aggregate'
TIMESTAMPS = ['intime', 'outtime', 'admittime', 'dischtime', 'deathtime', 'dod', 'corrected_death']


def classify_endpoint(frame):
    """Use exact intrahospital death ordering; retain unresolved dates as unknown."""
    exact = frame.deathtime.notna() & frame.deathtime.ge(frame.admittime)
    death = exact & frame.deathtime.ge(frame.intime) & frame.deathtime.le(frame.outtime)
    impossible = exact & frame.deathtime.lt(frame.intime)
    uncertain = (~exact & frame.dod.notna()
                 & frame.dod.dt.normalize().le(frame.outtime.dt.normalize())) | impossible
    status = pd.Series(2, index=frame.index, dtype='Int64')
    status.loc[death] = 1
    status.loc[uncertain] = pd.NA
    end = frame.outtime.where(~death, frame.deathtime).where(~uncertain)
    return status, end, exact, uncertain


def endpoint_tests():
    cases = [
        ('Death inside ICU', '2020-01-03 12:00', None, 1),
        ('Death at ICU exit', '2020-01-04 12:00', None, 1),
        ('Exact death after ICU exit', '2020-01-04 13:00', None, 2),
        ('No death recorded', None, None, 2),
        ('Date only after ICU exit date', None, '2020-01-05', 2),
        ('Date only on ICU exit date', None, '2020-01-04', None),
        ('Date only before ICU exit date', None, '2020-01-03', None),
        ('Invalid exact time and uncertain date', '2019-12-31', '2020-01-04', None),
        ('Death before ICU admission', '2020-01-01 12:00', None, None),
    ]
    test = pd.DataFrame(cases, columns=['case', 'deathtime', 'dod', 'expected'])
    for field, value in [('admittime', '2020-01-01'), ('intime', '2020-01-02'),
                         ('outtime', '2020-01-04 12:00')]:
        test[field] = pd.Timestamp(value)
    for field in ['deathtime', 'dod']:
        test[field] = pd.to_datetime(test[field], format='mixed')
    status, end, _, _ = classify_endpoint(test)
    for i, row in test.iterrows():
        expected = row.expected
        assert (pd.isna(status.loc[i]) if pd.isna(expected) else status.loc[i] == expected), row['case']
        assert pd.isna(end.loc[i]) == pd.isna(expected), row['case']
    # A death at a landmark has no subsequent follow-up and is excluded.
    assert not (pd.Timestamp('2020-01-03') > test.intime.iloc[0] + pd.Timedelta(days=1))
    return [dict(check=row['case'], status='PASS') for _, row in test.iterrows()]


def publish(frame, path):
    stage = path.with_suffix(path.suffix + '.tmp')
    frame.to_csv(stage, index=False)
    stage.replace(path)


def main():
    PRIVATE.mkdir(exist_ok=True)
    AGGREGATE.mkdir(exist_ok=True)
    tests = endpoint_tests()
    labels = pd.read_csv(SOURCE / 'data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    labels = labels.loc[labels.selected_k.eq(3), ['dataset', 'icu_stay_id', 'gbmtm_class', 'max_post_prob']]
    scales = pd.read_csv(SOURCE / 'aggregate/legacy_scale_parameters.csv').set_index('exposure')
    module_spec = importlib.util.spec_from_file_location('current_flows', WORK / 'display/prepare_flows.py')
    flows = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(flows)
    treatment = pd.read_csv(SOURCE / 'data_private/current_day1_treatments.csv')
    frames, populations, sources, descriptive, paths = [], [], [], [], []
    for db in ['mimiciv', 'mimiciii']:
        daily = pd.read_csv(SOURCE / f'data_private/all_icu_d1d7_daily_{db}.csv', low_memory=False)
        day1 = daily.loc[daily.day.eq(1)].copy()
        raw = pd.read_csv(SOURCE / f'data_private/{db}_raw_endpoints.csv')
        d = day1.merge(raw, on=['icu_stay_id', 'subject_id', 'encounter_id'], validate='one_to_one')
        assert len(d) == len(day1) == len(raw)
        d = d.merge(labels.loc[labels.dataset.eq(db)].drop(columns='dataset'), on='icu_stay_id',
                    how='left', validate='one_to_one')
        for field in TIMESTAMPS:
            d[field] = pd.to_datetime(d[field], format='mixed')
        status, end, exact, uncertain = classify_endpoint(d)
        d['icu_status'] = status
        d['icu_endpoint_time_days'] = (end - d.intime).dt.total_seconds() / 86400
        d['icu_exit_time_days'] = (d.outtime - d.intime).dt.total_seconds() / 86400
        d['icu_outcome'] = status.map({1: 'ICU death', 2: 'Live ICU discharge'}).fillna('Unknown ICU outcome')
        d['static_previous_eligible'] = d.followup_90d_days.gt(1)
        d['static_eligible'] = d.static_previous_eligible & d.icu_exit_time_days.gt(1)
        d['trajectory_previous_eligible'] = (d.gbmtm_class.notna() & d.icu_los_days.ge(7)
                                              & d.mortality_7d.eq(0) & d.followup_90d_days.gt(7))
        d['trajectory_eligible'] = d.trajectory_previous_eligible & d.icu_exit_time_days.gt(7)
        d['treatment_eligible'] = d.gbmtm_class.notna() & d.static_eligible
        for index in ['easix', 'siri']:
            value = d[f'{index}_day'].where(d[f'{index}_day'].gt(0))
            logged = np.log2(value)
            scale = scales.loc[index]
            d[f'log2_{index}_z'] = (logged.clip(scale.lower, scale.upper) - scale['mean']) / scale['sd']
        race = d.race_ethnicity.fillna('').str.upper().str.strip()
        d['race_group'] = np.select([race.str.contains('WHITE'), race.str.contains('BLACK|AFRICAN'),
            race.str.contains('HISPANIC|LATINO'), race.str.contains('ASIAN'),
            race.eq('') | race.str.contains('UNKNOWN|UNABLE|DECLINED')],
            ['White', 'Black', 'Hispanic', 'Asian', 'Unknown'], default='Other')
        d['sex_group'] = d.sex.str.upper().map({'M': 'Male', 'F': 'Female'}).fillna('Other')
        d['male'] = d.sex_group.eq('Male').astype(int)
        d['age_ge65'] = d.age.ge(65).astype(int)
        d['aki_day1'] = d.aki_stage_day1.ge(1).astype(float).where(d.aki_stage_day1.notna())
        for population, flag, previous in [
            ('static_day1_24h_landmark', 'static_eligible', 'static_previous_eligible'),
            ('day7_landmark', 'trajectory_eligible', 'trajectory_previous_eligible'),
            ('assigned_day1_treatment', 'treatment_eligible', 'treatment_eligible')]:
            x = d.loc[d[flag]]
            origin = 7 if population == 'day7_landmark' else 1
            known = x.icu_status.notna()
            assert x.loc[known, 'icu_endpoint_time_days'].gt(origin).all()
            populations.append(dict(dataset=db, analysis=population, origin_days=origin,
                previous_risk_n=int(d[previous].sum()), no_positive_icu_followup=int((d[previous] & ~d[flag]).sum()),
                n_landmark=len(x), unknown_icu_outcome=int((~known).sum()), n_endpoint_known=int(known.sum()),
                icu_deaths=int(x.icu_status.eq(1).sum()), live_discharges=int(x.icu_status.eq(2).sum())))
        sources.append(dict(dataset=db, parent_n=len(d), valid_exact_death=int(exact.sum()),
            invalid_exact_before_admission=int((d.deathtime.notna() & d.deathtime.lt(d.admittime)).sum()),
            icu_deaths=int(status.eq(1).sum()), live_discharges=int(status.eq(2).sum()),
            unknown_icu_outcome=int(uncertain.sum()), missing_outtime=int(d.outtime.isna().sum()),
            invalid_icu_interval=int(d.outtime.le(d.intime).sum())))
        for population, mask in [('parent', pd.Series(True, index=d.index)),
            ('fixed_k3_assigned', d.gbmtm_class.notna()), ('day7_landmark', d.trajectory_eligible)]:
            x = d.loc[mask]
            for cls in ['All', 1, 2, 3]:
                g = x if cls == 'All' else x.loc[x.gbmtm_class.eq(cls)]
                if not len(g):
                    continue
                deaths = int(g.icu_status.eq(1).sum())
                known_n = int(g.icu_status.notna().sum())
                descriptive.append(dict(dataset=db, population=population, trajectory_class=cls, n_total=len(g),
                    n_known=known_n, icu_deaths=deaths, live_discharges=int(g.icu_status.eq(2).sum()),
                    unknown=int(g.icu_status.isna().sum()), icu_mortality_pct_known=100 * deaths / known_n))
        assigned = d.loc[d.gbmtm_class.notna()].merge(treatment.loc[treatment.dataset.eq(db)].drop(columns='dataset'),
                on='icu_stay_id', validate='one_to_one')
        assigned['treatment_intensity'], assigned['fluid_cutoff_basis'] = flows.assign_intensity(assigned)
        n = len(assigned)
        for (cls, intensity), g in assigned.groupby(['gbmtm_class', 'treatment_intensity']):
            paths.append(dict(dataset=db, trajectory_class=int(cls), treatment_intensity=intensity,
                n=len(g), cohort_n=n, cohort_pct=100 * len(g) / n,
                icu_deaths=int(g.icu_status.eq(1).sum()), live_discharges=int(g.icu_status.eq(2).sum()),
                unknown_icu_outcome=int(g.icu_status.isna().sum())))
        treatment_columns = ['icu_stay_id', 'oasis_without_mv', 'nee_any_24h', 'nee_max_24h',
            'input_documented', 'output_documented', 'fluid_balance_ml_24h', 'fluid_balance_ml_kg_24h',
            'weight_kg', 'documented_mv_24h', 'documented_crrt_24h', 'treatment_intensity']
        d = d.merge(assigned[treatment_columns], on='icu_stay_id', how='left', validate='one_to_one')
        for flag, reference in [('static_previous_eligible', 'static_alive_24h'),
                                ('trajectory_previous_eligible', 'trajectory_alive_icu_resident_day7')]:
            previous = pd.read_csv(SOURCE / 'aggregate/landmark_population_flow.csv')
            expected = int(previous.loc[previous.dataset.eq(db) & previous.stage.eq(reference), 'n'].iloc[0])
            assert d[flag].sum() == expected, (db, flag)
        frames.append(d)
    all_rows = pd.concat(frames, ignore_index=True)
    assert not all_rows.duplicated(['dataset', 'icu_stay_id']).any()
    publish(all_rows, PRIVATE / 'icu_mortality_analysis.csv')
    publish(pd.DataFrame(populations), AGGREGATE / 'population_flow.csv')
    publish(pd.DataFrame(sources), AGGREGATE / 'endpoint_source_audit.csv')
    publish(pd.DataFrame(descriptive), AGGREGATE / 'outcome_counts.csv')
    path_table = pd.DataFrame(paths)
    for column in ['icu_deaths', 'live_discharges', 'unknown_icu_outcome']:
        path_table[column + '_cohort_pct'] = 100 * path_table[column] / path_table.cohort_n
    assert path_table.n.eq(path_table.icu_deaths + path_table.live_discharges + path_table.unknown_icu_outcome).all()
    publish(path_table, AGGREGATE / 'class_treatment_outcome_counts.csv')
    record = dict(created=datetime.now(timezone.utc).isoformat(), status='PASS',
        endpoint_definition_tests=tests, patient_rows=len(all_rows),
        article_outputs_updated=False, classifier_or_trajectory_refitted=False,
        sources='Current saved endpoint and fixed-class inputs; no new database extraction')
    (AGGREGATE / 'preparation_checks.json').write_text(json.dumps(record, indent=2) + '\n')
    print(pd.DataFrame(populations).to_string(index=False))


if __name__ == '__main__':
    main()
