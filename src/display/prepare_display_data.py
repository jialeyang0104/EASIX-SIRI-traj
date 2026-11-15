#!/usr/bin/env python3
"""Publish aggregate display data from fixed labels and existing local caches.

Uses current fixed-projection labels and strict 0–24 h multisystem means.
"""
from pathlib import Path
import argparse
import json
import re

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
W = HERE.parent
P = W.parent
A = W / 'endpoints/aggregate'
S = HERE.parents[1] / 'data/aggregate/display'
FIXED = P/'results/trajectory/tables'
PROFILE = P/'data/private/profile_tables'
COLORS = {'1': '#4C78A8', '2': '#59A14F', '3': '#E15759'}
DB = {'mimiciv': 'MIMIC-IV', 'mimiciii': 'MIMIC-III'}
VARIABLES = [
    ('age', 'Age', 'years', 'Demography', False),
    ('oasis_day1', 'OASIS', 'points', 'Severity', False),
    ('sofa_day1', 'SOFA', 'points', 'Severity', False),
    ('apsiii_day1', 'APS III', 'points', 'Severity', False),
    ('easix_day', 'EASIX', 'index', 'Endothelial stress', True),
    ('ldh_u_l', 'LDH', 'U/L', 'Endothelial stress', True),
    ('creatinine_mg_dl', 'Creatinine', 'mg/dL', 'Renal', True),
    ('platelets_109_l', 'Platelets', '10⁹/L', 'Haematology', True),
    ('siri_day', 'SIRI', 'index', 'Inflammation', True),
    ('wbc_109_l', 'White blood cells', '10⁹/L', 'Inflammation', True),
    ('neutrophils_harmonised', 'Neutrophils', '10⁹/L', 'Inflammation', True),
    ('monocytes_harmonised', 'Monocytes', '10⁹/L', 'Inflammation', True),
    ('lymphocytes_harmonised', 'Lymphocytes', '10⁹/L', 'Inflammation', True),
]

def publish(frame, name):
    path = S / (name + '.csv')
    stage = path.with_suffix('.csv.tmp')
    frame.to_csv(stage, index=False)
    stage.replace(path)
    return pd.read_csv(path)

def day1(db):
    path = W / 'endpoints/data_private' / f'all_icu_d1d7_daily_{db}.csv'
    d = pd.read_csv(path, low_memory=False)
    d = d.loc[d.day.eq(1)].copy()
    assert d.icu_stay_id.is_unique
    # Counts have the same units in both databases. Reconstructed components use
    # the same daily WBC and differential summaries as the source SIRI contract.
    use_absolute = d[['neutrophils_abs', 'monocytes_abs', 'lymphocytes_abs']].notna().all(axis=1)
    if db == 'mimiciii':
        use_absolute[:] = False
    for component in ['neutrophils', 'monocytes', 'lymphocytes']:
        absolute = d[component + '_abs']
        reconstructed = d.wbc_109_l * d[component + '_pct'] / 100
        d[component + '_harmonised'] = absolute.where(use_absolute, reconstructed)
    d['dataset'] = db
    d['male']=d.sex.eq('M').astype(float).where(d.sex.notna())
    means = pd.read_csv(W / 'endpoints/data_private' / f'{db}_day1_system_means.csv')
    return d.merge(means, on='icu_stay_id', validate='one_to_one')

MANAGEMENT_METRICS = [
    ('nee_any','Any recorded NEE','proportion',True),
    ('nee_high','Maximum NEE >0.3 µg/kg/min','proportion',True),
    ('fluid_kg_high','Net fluid balance >30 mL/kg','proportion',True),
    ('mv','Recorded invasive mechanical ventilation','proportion',True),
    ('crrt','Recorded CRRT','proportion',True),
    ('map','MAP · mean','mmHg',False),('urine','Urine output','mL/kg/h',False),
    ('lactate','Lactate · mean','mmol/L',False)]

def treatment_values(d):
    d=d.copy()
    d['nee_any']=d.nee_any_24h
    d['nee_high']=d.nee_max_24h.gt(.3).astype(float).where(d.nee_max_24h.notna())
    valid=d.input_documented.eq(1)&d.output_documented.eq(1)&d.weight_kg.between(20,300)
    d['fluid_kg_high']=d.fluid_balance_ml_kg_24h.gt(30).astype(float).where(valid)
    d['mv']=d.documented_mv_24h;d['crrt']=d.documented_crrt_24h
    return d

def display_p(groups,binary=False):
    if not all(len(x) for x in groups):return np.nan
    if binary:
        cells=np.array([[x.sum(),len(x)-x.sum()] for x in groups]);cells=cells[:,cells.sum(axis=0)>0]
        return stats.chi2_contingency(cells)[1] if cells.shape[1]>1 else 1.0
    if len(set(np.concatenate(groups)))<2:return 1.0
    return stats.kruskal(*groups).pvalue

def baseline(cohorts):
    """Publication parent/class tables with explicit observed denominators."""
    labels=pd.read_csv(W/'endpoints/data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    tx=treatment_values(pd.read_csv(W/'endpoints/data_private/current_parent_day1_treatments.csv',low_memory=False))
    data={db:d.merge(tx[tx.dataset.eq(db)][['icu_stay_id']+[m[0] for m in MANAGEMENT_METRICS[:5]]],on='icu_stay_id',validate='one_to_one') for db,d in cohorts.items()}
    icu=pd.read_csv(W/'icu/data_private/icu_mortality_analysis.csv',usecols=['dataset','icu_stay_id','icu_status'])
    for db,d in data.items():
        d=d.merge(icu.loc[icu.dataset.eq(db),['icu_stay_id','icu_status']],on='icu_stay_id',validate='one_to_one')
        d['icu_mortality']=d.icu_status.eq(1).astype(float).where(d.icu_status.notna())
        d['icu_outcome_unknown']=d.icu_status.isna().astype(float)
        data[db]=d
    mapping=pd.read_csv(S/'display_variables.csv')
    specs=[('age','Age, years',False),('male','Male sex',True),('white','White race/ethnicity',True),('sepsis','Sepsis',True),('oasis_day1','OASIS',False),('sofa_day1','SOFA',False),('apsiii_day1','APS III',False),('sirs_day1','SIRS score',False),('mechanical_ventilation_day1','Mechanical ventilation, Day 1',True),('vasopressor_day1','Vasopressor, Day 1',True),('crrt_day1','CRRT, Day 1',True),('aki','AKI, Day 1',True)]
    specs += [('ldh_u_l','LDH · maximum, U/L',False),('creatinine_mg_dl','Creatinine · maximum, mg/dL',False),('platelets_109_l','Platelets · minimum, 10⁹/L',False),('wbc_109_l','White blood cells · maximum, 10⁹/L',False),('neutrophils_harmonised','Neutrophils, 10⁹/L',False),('monocytes_harmonised','Monocytes, 10⁹/L',False),('lymphocytes_harmonised','Lymphocytes, 10⁹/L',False),('neutrophils_pct','Neutrophils, %',False),('monocytes_pct','Monocytes, %',False),('lymphocytes_pct','Lymphocytes, %',False)]
    specs += [(r.variable,r.label+(' · mean' if r.variable.endswith('_mean') else ' · daily index')+', '+r.unit,False) for r in mapping.itertuples()]
    specs += [(m,l,True) for m,l,u,b in MANAGEMENT_METRICS[:5]]
    specs += [('icu_los_days','ICU length of stay, days',False),('hospital_los_days','Hospital length of stay, days',False),('mortality_28d','28-day mortality',True),('icu_mortality','ICU mortality (known outcomes)',True),('icu_outcome_unknown','ICU outcome unknown',True)]
    for db,d in data.items():
        d['white']=d.race_ethnicity.str.upper().str.startswith('WHITE').astype(float).where(d.race_ethnicity.notna());d['aki']=d.aki_stage_day1.gt(0).astype(float).where(d.aki_stage_day1.notna())
    def cell(x,binary,variable):
        x=x.replace([np.inf,-np.inf],np.nan).dropna();n=len(x)
        if not n:return 'NA [observed n=0]'
        if binary:return f'{int(x.sum())}/{n} ({100*x.mean():.1f}%)'
        digits=2 if variable=='ph_mean' else 1
        return f'{x.median():.{digits}f} ({x.quantile(.25):.{digits}f}–{x.quantile(.75):.{digits}f}) [observed n={n}]'
    def pv(p):return '' if pd.isna(p) else '<0.001' if p<.001 else f'{p:.3f}'
    rows=[dict(Characteristic='Patients',Type='count',**{DB[db]:str(len(d)) for db,d in data.items()},**{'P value':''})];classes=[];counts=[]
    seen=set()
    for var,label,b in specs:
        if var in seen:continue
        seen.add(var);xs=[d[var].dropna() for d in data.values()]
        rows.append(dict(Characteristic=label,Type='binary' if b else 'continuous',**{DB[db]:cell(d[var],b,var) for db,d in data.items()},**{'P value':pv(display_p(xs,b))}))
        for db,d in data.items():
            z=d.merge(labels[labels.dataset.eq(db)][['icu_stay_id','gbmtm_class']],on='icu_stay_id',validate='one_to_one');groups=[z.loc[z.gbmtm_class.eq(c),var].dropna() for c in [1,2,3]]
            classes.append(dict(Dataset=DB[db],Characteristic=label,**{f'Class {c}':cell(x,b,var) for c,x in zip([1,2,3],groups)},**{'P value':pv(display_p(groups,b))}))
            for pop,a in [('Parent cohort',d),('Trajectory cohort',z)]:
                counts.append(dict(dataset=db,population=pop,variable=var,label=label,stays=len(a),observed=a[var].notna().sum(),missing=a[var].isna().sum(),missing_pct=100*a[var].isna().mean()))
    for db,d in data.items():
        z=d.merge(labels[labels.dataset.eq(db)][['icu_stay_id','gbmtm_class']],on='icu_stay_id');classes.insert(0,dict(Dataset=DB[db],Characteristic='Patients',**{f'Class {c}':str(z.gbmtm_class.eq(c).sum()) for c in [1,2,3]},**{'P value':''}))
    publish(pd.DataFrame(rows)[['Characteristic','Type','MIMIC-IV','MIMIC-III','P value']],'publication_parent_baseline')
    publish(pd.DataFrame(classes)[['Dataset','Characteristic','Class 1','Class 2','Class 3','P value']],'publication_class_baseline')
    publish(pd.DataFrame(rows),'table_1_expanded_baseline')
    publish(pd.DataFrame(counts),'baseline_availability')

def phenotypes(cohorts):
    labels = pd.read_csv(W / 'endpoints/data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    ds = {}
    for db, d in cohorts.items():
        lab = labels.loc[labels.dataset.eq(db), ['icu_stay_id', 'gbmtm_class']]
        ds[db] = d.merge(lab, on='icu_stay_id', validate='one_to_one')
    reference = ds['mimiciv']
    transformed, parameters = {}, []
    for variable, label, unit, system, log in VARIABLES:
        ref = reference[variable].astype(float)
        if log:
            ref = np.log2(ref.where(ref.gt(0)))
        low, high = ref.quantile([.01, .99])
        clipped = ref.clip(low, high)
        mean, sd = clipped.mean(), clipped.std(ddof=1)
        parameters.append(dict(variable=variable, label=label, unit=unit, system=system,
                               log2=log, lower=low, upper=high, mean=mean, sd=sd,
                               reference_n=ref.notna().sum(), reference_population='MIMIC-IV fixed trajectory cohort Day 1'))
        for db, d in ds.items():
            raw = d[variable].astype(float)
            value = np.log2(raw.where(raw.gt(0))) if log else raw
            transformed[db, variable] = (value.clip(low, high)-mean)/sd
    publish(pd.DataFrame(parameters), 'phenotype_scaling')
    baseline_rows, radar_rows, association_rows = [], [], []
    for db, d in ds.items():
        for variable, label, unit, system, log in VARIABLES:
            z = transformed[db, variable]
            ref = np.sort(reference[variable].dropna().to_numpy())
            raw = d[variable]
            percentiles = pd.Series(100*(np.searchsorted(ref, raw, side='left') +
                                          np.searchsorted(ref, raw, side='right'))/(2*len(ref)), index=d.index)
            percentiles[raw.isna()] = np.nan
            for c in [1, 2, 3]:
                mask = d.gbmtm_class.eq(c)
                x = z.loc[mask].dropna()
                r = raw.loc[mask].dropna()
                n = len(x)
                ci = stats.t.ppf(.975, n-1)*x.std(ddof=1)/np.sqrt(n)
                base = dict(scope=db, state=c, variable=variable, label=label, system=system,
                            stays=int(mask.sum()), observed=n, missing=int(mask.sum())-n,
                            unit=unit, source='projected labels + Day 1 cache',
                            result_batch='fixed K3 d2 pruning; current landmark endpoints')
                baseline_rows.append(dict(**base, raw_mean=r.mean(), median=r.median(),
                                          q25=r.quantile(.25), q75=r.quantile(.75),
                                          mean_z=x.mean(), ci_low=x.mean()-ci, ci_high=x.mean()+ci,
                                          interval_method='t mean CI conditional on fixed labels'))
                radar_rows.append(dict(**base, mean_percentile=percentiles.loc[mask].mean(), reference_n=len(ref)))
                valid = z.notna()
                rr, pp = stats.pointbiserialr(mask.loc[valid].astype(int), z.loc[valid])
                association_rows.append(dict(**base, paired_n=int(valid.sum()), r=rr,
                                             status='estimable', contrast='class vs all other classes'))
    publish(pd.DataFrame(baseline_rows), 'baseline_comparison')
    publish(pd.DataFrame(radar_rows), 'subtype_radar')
    publish(pd.DataFrame(association_rows), 'system_associations')

def system_radar():
    """Give each clinical system one equally weighted axis from observed means."""
    radar = pd.read_csv(S/'subtype_radar.csv')
    mapping = pd.read_csv(S/'display_variables.csv')
    mapping = mapping.loc[mapping.system.ne('Index components')].copy()
    mapping['variables_in_system'] = mapping.groupby('system').variable.transform('size')
    mapping['feature_weight'] = 1 / mapping.variables_in_system
    mapping['source'] = 'subtype_radar.csv: observed-variable mean reference percentiles'
    publish(mapping, 'system_radar_variables')
    rows = []
    for (db, c, system), group in radar.loc[radar.variable.isin(mapping.variable)].groupby(
            ['scope', 'state', 'system'], sort=False):
        expected = mapping.loc[mapping.system.eq(system), 'variable'].tolist()
        assert set(group.variable) == set(expected) and group.variable.is_unique
        assert group.mean_percentile.notna().all()
        assert group.stays.nunique() == 1
        rows.append(dict(scope=db, state=c, system=system,
                         stays=int(group.stays.iloc[0]), variables_in_system=len(group),
                         mean_percentile=group.mean_percentile.mean(),
                         minimum_variable_n=int(group.observed.min()),
                         maximum_variable_n=int(group.observed.max()),
                         variables=';'.join(expected),
                         aggregation='Equal mean of variable-specific observed mean reference percentiles',
                         reference='MIMIC-IV trajectory-cohort observed-value midranks',
                         source='subtype_radar.csv; system_radar_variables.csv'))
    result = pd.DataFrame(rows)
    assert len(result) == 54 and result.system.nunique() == 9
    publish(result, 'subtype_system_radar')

def management():
    tx=treatment_values(pd.read_csv(W/'endpoints/data_private/current_parent_day1_treatments.csv',low_memory=False))
    labels=pd.read_csv(W/'endpoints/data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    rows=[]
    for db in DB:
        d=day1(db).merge(tx[tx.dataset.eq(db)].drop(columns='dataset'),on='icu_stay_id',validate='one_to_one').merge(labels[labels.dataset.eq(db)][['icu_stay_id','gbmtm_class']],on='icu_stay_id',validate='one_to_one')
        d['map']=d.map_mean;d['lactate']=d.lactate_mean
        # Strict 0–24 h urine volume from the four six-hour intervals; denominator 24 h.
        prof=pd.read_csv(W/f'endpoints/data_private/{db}_current_six_hour_management.csv')
        urine=prof.groupby('icu_stay_id').urine_output_ml.sum(min_count=1)
        d['urine']=d.icu_stay_id.map(urine)/24/d.weight_kg
        for metric,label,unit,b in MANAGEMENT_METRICS:
            groups=[d.loc[d.gbmtm_class.eq(c),metric].dropna() for c in [1,2,3]];p=display_p(groups,b)
            for c,x in zip([1,2,3],groups):
                center,lo,hi=management_interval(x,b)
                rows.append(dict(dataset=db,trajectory_class=c,metric=metric,metric_label=label,unit=unit,binary=b,center=center,low=lo,high=hi,n_observed=len(x),cohort_n=int(d.gbmtm_class.eq(c).sum()),p_value=p,interval_method='Wilson 95% CI' if b else 'median and IQR',source='Strict 0–24 h system cache arithmetic mean (MAP/lactate); sum urine volume /24 h /valid weight; current recorded treatments'))
    publish(pd.DataFrame(rows),'management_day1_summary')

def management_interval(x,binary):
    n=len(x)
    if not n:return np.nan,np.nan,np.nan
    center=x.mean() if binary else x.median()
    if not binary:return center,x.quantile(.25),x.quantile(.75)
    z=stats.norm.ppf(.975);den=1+z*z/n;mid=(center+z*z/(2*n))/den;half=z*np.sqrt(center*(1-center)/n+z*z/(4*n*n))/den
    return center,mid-half,mid+half

def supplementary_management():
    rows=[]
    for db in DB:
        d=pd.read_csv(W/f'endpoints/data_private/{db}_current_six_hour_management.csv').sort_values(['icu_stay_id','interval_index'])
        d['nee_high']=d.nee_max.gt(.3).astype(float).where(d.nee_max.notna())
        g=d.groupby("icu_stay_id")
        # Cumulative recorded volumes remain usable after an empty later interval.
        inp=d.input_ml.fillna(0).groupby(d.icu_stay_id).cumsum();out=d.output_ml.fillna(0).groupby(d.icu_stay_id).cumsum()
        valid=g.input_events.cumsum().gt(0)&g.output_events.cumsum().gt(0)&d.weight_kg.between(20,300)
        d['fluid_kg_high']=((inp-out)/d.weight_kg).gt(30).astype(float).where(valid)
        d['map']=d.map_median;d['lactate']=d.lactate_median;d['urine']=d.urine_output_ml_kg_h
        for (c,k),group in d.groupby(['trajectory_class','interval_index']):
            for metric,label,unit,b in MANAGEMENT_METRICS:
                x=group[metric].dropna();center,lo,hi=management_interval(x,b)
                rows.append(dict(dataset=db,trajectory_class=int(c),trajectory_label=f'Class {c}',metric=metric,metric_label=label.replace(' · mean',' · interval median'),unit=unit,binary=b,time_start_h=6*(k-1),time_end_h=6*k,time_label=f'{6*(k-1)}–{6*k} h',center=center,p25=np.nan if b else lo,p75=np.nan if b else hi,ci_lower=lo if b else np.nan,ci_upper=hi if b else np.nan,n_observed=len(x),n_at_risk=len(group),cohort_n=len(group),missing_pct=100*(1-len(x)/len(group)),center_statistic='proportion' if b else 'median',spread_type='Wilson 95% CI' if b else 'IQR',source_granularity='Fixed classes; interval NEE maximum and recorded MV/CRRT; cumulative admission-to-endpoint fluid balance >30 mL/kg with recorded input/output and valid weight; interval physiology'))
    result=pd.DataFrame(rows);assert len(result)==192 and result.groupby(['dataset','metric','trajectory_class']).size().eq(4).all()
    publish(result,'supp_management_profiles')

def flow():
    d = pd.read_csv(A/'all_icu_flow_diagram_counts.csv')
    exclusions = {
        'First ICU stay per patient': 'Subsequent ICU stays',
        'Age >=18 years': 'Age <18 years',
        'ICU length of stay >=24 h': 'ICU stay <24 h',
        'Hospital length of stay >=24 h': 'Hospital stay <24 h',
    }
    config = dict(source='src/endpoints/aggregate/all_icu_flow_diagram_counts.csv',
                  groups=[dict(label='Development cohort', cohorts=['mimiciv'], palette='primary'),
                          dict(label='Historical cross-version cohort', cohorts=['mimiciii'], palette='external')], cohorts=[],
                  style=dict(width_mm=183, dpi=300, font_size_pt=8, exclusion_font_size_pt=8))
    static = pd.read_csv(A/'all_icu_cox_static_effects.csv')
    pop = pd.read_csv(A/'landmark_population_flow.csv')
    analysis_rows=[]
    for db in DB:
        seq = d.loc[d.dataset.eq(db) & d.stage_type.eq('sequential')].sort_values('step_order')
        nodes=[dict(label=r.step.replace('>=','≥'),n=int(r.n),unit='stays' if i==0 else 'patients') for i,r in enumerate(seq.itertuples())]
        nodes[-1]['label']='Eligible parent cohort'
        ex=[dict(label=exclusions[r.step],n=int(r.excluded_from_previous)) for r in list(seq.itertuples())[1:]]
        config['cohorts'].append(dict(id=db,label=DB[db],nodes=nodes,exclusions=ex))
        available=d[(d.dataset==db)&(d.stage_type=='branch')].set_index('stage_id')['n']
        current=static[(static.dataset==db)&(static.analysis=='static_day1_24h_landmark')]
        assert current.groupby('model').n.nunique().eq(1).all(), 'Static sample size differs across horizons'
        models=current.groupby('model')['n'].first().to_dict()
        assigned=int(pop[(pop.dataset==db)&(pop.stage=='fixed_k3_assigned')].n.iloc[0]); day7=int(pop[(pop.dataset==db)&(pop.stage=='trajectory_alive_icu_resident_day7')].n.iloc[0])
        config['cohorts'][-1]['analysis']=dict(static_available={k:int(available['static_'+k]) for k in ['easix','siri','joint']},static_models={k:int(models['static_'+k]) for k in ['easix','siri','joint']},trajectory_assigned=assigned,trajectory_day7=day7)
        analysis_rows += [dict(dataset=db,population='static_available_easix',n=int(available['static_easix']),source='all_icu_flow_diagram_counts.csv'),dict(dataset=db,population='static_available_siri',n=int(available['static_siri']),source='all_icu_flow_diagram_counts.csv'),dict(dataset=db,population='static_available_joint',n=int(available['static_joint']),source='all_icu_flow_diagram_counts.csv'),dict(dataset=db,population='static_model_easix',n=int(models['static_easix']),source='all_icu_cox_static_effects.csv'),dict(dataset=db,population='static_model_siri',n=int(models['static_siri']),source='all_icu_cox_static_effects.csv'),dict(dataset=db,population='static_model_joint',n=int(models['static_joint']),source='all_icu_cox_static_effects.csv'),dict(dataset=db,population='trajectory_assigned',n=assigned,source='landmark_population_flow.csv'),dict(dataset=db,population='trajectory_day7',n=day7,source='landmark_population_flow.csv')]
    (S/'figure_1_cohort_flow.json').write_text(json.dumps(config,ensure_ascii=False,indent=2)+'\n')
    publish(pd.DataFrame(analysis_rows),'figure_1_analysis_populations')
    publish(pop,'landmark_populations')

def trajectory_comparisons(cohorts):
    labels=pd.read_csv(W/'endpoints/data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    rows=[]
    for db,d in cohorts.items():
        daily=pd.read_csv(W/'endpoints/data_private'/f'all_icu_d1d7_daily_{db}.csv',low_memory=False)
        x=daily.merge(labels.loc[labels.dataset.eq(db),['icu_stay_id','gbmtm_class']],on='icu_stay_id',validate='many_to_one')
        for index,var in [('log2_easix_interp_z','easix_day'),('log2_siri_interp_z','siri_day')]:
            pvals=[]; pending=[]
            for day in range(1,8):
                groups=[x.loc[x.gbmtm_class.eq(c)&x.day.eq(day),var].dropna() for c in [1,2,3]]
                result=stats.kruskal(*groups);pvals.append(result.pvalue);pending.append((day,result,groups))
            q=np.empty(7); order=np.argsort(pvals); running=1
            for pos in range(6,-1,-1):
                i=order[pos]; running=min(running,pvals[i]*7/(pos+1));q[i]=running
            for (day,result,groups),pv,qv in zip(pending,pvals,q):
                stars='***' if qv<.001 else '**' if qv<.01 else '*' if qv<.05 else 'ns'
                rows.append(dict(dataset=db,index=index,day=day,statistic=result.statistic,p_value=pv,p_value_bh=qv,n=sum(len(g) for g in groups),n_class1=len(groups[0]),n_class2=len(groups[1]),n_class3=len(groups[2]),stars=stars,test='Kruskal-Wallis (2 df)',interpretation='descriptive comparison of observed values conditional on fixed labels'))
    result=pd.DataFrame(rows)
    assert len(result)==28
    publish(result,'trajectory_daily_comparisons')
    publish(result.loc[result.day.eq(1)],'trajectory_day1_comparisons')

def configure_systems(cohorts):
    global VARIABLES
    metadata = pd.read_csv(A/'day1_system_variable_dictionary.csv').drop_duplicates('variable')
    labels = pd.read_csv(W/'endpoints/data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    included=[];counts=[]
    # Include a common set of measures with <70% missing in both labelled cohorts.
    # Index-related measurements are exempt from this threshold.
    for r in metadata.itertuples():
        rates=[]
        for db,d in cohorts.items():
            labelled=d[d.icu_stay_id.isin(labels.loc[labels.dataset.eq(db),'icu_stay_id'])]
            rates.append(labelled[r.variable].isna().mean())
        if r.index_component or max(rates)<.70:
            included.append((r.variable,r.label,r.unit,r.system,r.variable not in ['ph_mean','baseexcess_mean','temperature_mean','gcs_mean','spo2_mean']))
    # Keep both index definitions alongside the independent system means.
    VARIABLES=[('easix_day','EASIX','index','Index components',True),('siri_day','SIRI','index','Index components',True)]+included
    for db,d in cohorts.items():
        labelled=d[d.icu_stay_id.isin(labels.loc[labels.dataset.eq(db),'icu_stay_id'])]
        specs=[(r.variable,r.label,r.unit,r.system) for r in metadata.itertuples()]
        specs += [(v,l,u,sy) for v,l,u,sy,_ in ORIGINAL_VARIABLES]+[
            ('sex','Sex','category','Demography'),('male','Male sex','binary','Demography'),
            ('race_ethnicity','Race/ethnicity','category','Demography'),('sepsis','Sepsis','binary','Clinical context'),
            ('aki_stage_day1','AKI stage','stage','Organ support'),
            ('sirs_day1','SIRS','points','Severity'),('mechanical_ventilation_day1','Mechanical ventilation','binary','Organ support'),
            ('vasopressor_day1','Any vasopressor','binary','Organ support'),('crrt_day1','CRRT','binary','Organ support'),
            ('neutrophils_pct','Neutrophils','%','Index components'),('monocytes_pct','Monocytes','%','Index components'),('lymphocytes_pct','Lymphocytes','%','Index components')]
        seen=set()
        for v,l,u,sy in specs:
            if v in seen:continue
            seen.add(v)
            for population,x in [('Parent cohort',d),('Trajectory cohort',labelled)]:
                n=int(x[v].notna().sum())
                counts.append(dict(dataset=db,population=population,variable=v,label=l,unit=u,system=sy,stays=len(x),observed=n,missing=len(x)-n,missing_pct=100*(len(x)-n)/len(x),figure4_included=v in [a[0] for a in VARIABLES]))
    publish(pd.DataFrame(counts),'baseline_availability')
    publish(pd.DataFrame([dict(variable=v,label=l,unit=u,system=sy,log2=log) for v,l,u,sy,log in VARIABLES]),'display_variables')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--figure-sources-only', action='store_true')
    parser.add_argument('--system-radar-only', action='store_true',
                        help='Aggregate current variable percentiles without rebuilding other sources')
    parser.add_argument('--treatment-forest-only',action='store_true',
                        help='Bind current two-database cutoff treatment models without refitting')
    parser.add_argument('--supp-management-only',action='store_true',help='Add III to the retained six-hour management profile')
    args = parser.parse_args()
    if args.supp_management_only:
        supplementary_management()
        print('Published both databases for six-hour management profiles.')
        return
    if args.treatment_forest_only:
        treatment_forest()
        print('Published current 28-day/ICU cutoff treatment displays.')
        return
    if args.figure_sources_only:
        cohorts = {db: day1(db) for db in DB}
        flow()
        trajectory_comparisons(cohorts)
        print('Published cohort analysis populations and daily trajectory comparisons.')
        return
    if args.system_radar_only:
        system_radar()
        print('Published nine-system radar aggregates and equal-weight variable mapping.')
        return
    S.mkdir(exist_ok=True, parents=True)
    cohorts = {db: day1(db) for db in DB}
    configure_systems(cohorts)
    baseline(cohorts)
    phenotypes(cohorts)
    system_radar()
    management()
    supplementary_management()
    flow()
    # Legacy six-metric p-values superseded by management_day1_summary.csv.
    trajectory_comparisons(cohorts)
    treatment_forest()
    print('Published expanded baseline, strict system profiles, management P values and harmonised cohort flow.')

if __name__ == '__main__':
    main()
