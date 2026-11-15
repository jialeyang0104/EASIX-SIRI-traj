#!/usr/bin/env python3
"""Aggregate fixed-class ICU disposition and first-day treatment/outcome flows."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
W=HERE.parent;S=HERE.parents[1]/'data/aggregate/display'
PRIVATE=W/'endpoints/data_private'
from figure_support import validate_flow

def publish(frame,name):
    target=S/(name+'.csv');stage=target.with_suffix('.csv.tmp')
    frame.to_csv(stage,index=False);stage.replace(target)

def aggregate_paths(paths,meta):
    stages=sorted(paths)
    assert all(v.index.equals(paths[stages[0]].index) for v in paths.values())
    n=len(paths[stages[0]]);nodes=[];links=[]
    for stage in stages:
        for category,count in paths[stage].value_counts().items():
            percent=100*int(count)/n
            nodes.append(dict(**meta,hour=stage,category=category,stays=int(count),cohort_stays=n,
                panel_percent=percent,percentage_label=f'{percent:.1f}%'))
    origin=paths[stages[0]]
    for a,b in zip(stages[:-1],stages[1:]):
        frame=pd.DataFrame(dict(origin=origin,from_category=paths[a],to_category=paths[b]))
        for key,count in frame.groupby(['origin','from_category','to_category'],dropna=False).size().items():
            links.append(dict(**meta,from_hour=a,to_hour=b,origin=key[0],from_category=key[1],to_category=key[2],stays=int(count),cohort_stays=n))
    nn=pd.DataFrame(nodes);ll=pd.DataFrame(links)
    validate_flow(nn,ll,stages,categories=list(dict.fromkeys(nn.category)))
    return nn,ll

def disposition(labels):
    nodes=[];links=[];summaries=[]
    for db in ['mimiciv','mimiciii']:
        daily=pd.read_csv(PRIVATE/f'all_icu_d1d7_daily_{db}.csv',low_memory=False)
        daily=daily[daily.day.eq(1)]
        raw=pd.read_csv(PRIVATE/f'{db}_raw_endpoints.csv')
        d=labels[labels.dataset.eq(db)][['icu_stay_id','gbmtm_class']].merge(daily,on='icu_stay_id',validate='one_to_one').merge(raw,on=['icu_stay_id','subject_id','encounter_id'],validate='one_to_one')
        for field in ['intime','outtime','admittime','deathtime','dod','corrected_death']:
            d[field]=pd.to_datetime(d[field],format='mixed')
        exact=d.deathtime.notna()&d.deathtime.ge(d.admittime)
        icu_death=exact&d.deathtime.ge(d.intime)&d.deathtime.le(d.outtime)
        ambiguous_exit=~exact&d.dod.notna()&d.dod.dt.normalize().le(d.outtime.dt.normalize())
        paths={}
        for day in [1,3,5,7]:
            point=d.intime+pd.to_timedelta(day,unit='d')
            status=('Class '+d.gbmtm_class.astype(str)).copy()
            exited=d.outtime.lt(point)
            status[exited]='live_exit'
            status[exited&ambiguous_exit]='unknown'
            status[icu_death&d.deathtime.le(point)]='icu_death'
            # Date-only deaths before the checkpoint cannot be labelled resident/alive.
            status[~exited&~exact&d.corrected_death.le(point)]='unknown'
            paths[day]=status
            summaries.append(dict(dataset=db,day=day,cohort_n=len(d),
                resident=sum(status.eq(f'Class {c}').sum() for c in [1,2,3]),
                live_icu_exit=int(status.eq('live_exit').sum()),icu_death=int(status.eq('icu_death').sum()),
                unknown=int(status.eq('unknown').sum())))
        nn,ll=aggregate_paths(paths,dict(dataset=db,flow='icu_disposition',horizon=7))
        nodes.append(nn);links.append(ll)
        day7=next(r for r in summaries if r['dataset']==db and r['day']==7)
        expected=int(((d.icu_los_days>=7)&d.mortality_7d.eq(0)&d.followup_90d_days.gt(7)).sum())
        assert day7['resident']==expected,(db,day7,expected)
    publish(pd.concat(nodes,ignore_index=True),'supp_disposition_nodes')
    publish(pd.concat(links,ignore_index=True),'supp_disposition_links')
    publish(pd.DataFrame(summaries),'supp_disposition_summary')

def assign_intensity(d):
    """Organ-support-first hierarchy confirmed by the user; evaluate top down."""
    mv=d.documented_mv_24h.eq(1);crrt=d.documented_crrt_24h.eq(1)
    nee=d.nee_max_24h
    assert nee.notna().all(), 'Unquantifiable NEE requires source resolution before applying the confirmed hierarchy'
    io=d.input_documented.eq(1)&d.output_documented.eq(1)&np.isfinite(d.fluid_balance_ml_24h)
    weight=d.weight_kg.between(20,300)
    fluid_high=pd.Series(pd.NA,index=d.index,dtype='boolean')
    fluid_high.loc[io&weight]=d.loc[io&weight,'fluid_balance_ml_kg_24h'].gt(30)
    fluid_high.loc[io&~weight]=d.loc[io&~weight,'fluid_balance_ml_24h'].gt(2000)
    intensity=pd.Series('Intensity unknown',index=d.index)
    maximum=crrt|(mv&nee.gt(.3));intensity[maximum]='Max'
    high=~maximum&(nee.gt(.3)|(mv&nee.gt(0)&nee.le(.3)));intensity[high]='High'
    medium=~maximum&~high&((mv&nee.eq(0))|(~mv&nee.gt(0)&nee.le(.3))|fluid_high.fillna(False))
    intensity[medium]='Medium'
    low=~mv&~crrt&nee.eq(0)&fluid_high.eq(False).fillna(False);intensity[low]='Low'
    unknown=~mv&~crrt&nee.eq(0)&fluid_high.isna()
    assert intensity.eq('Intensity unknown').equals(unknown)
    basis=pd.Series('Missing fluid documentation',index=d.index)
    basis[io&weight]='30 mL/kg';basis[io&~weight]='2,000 mL (invalid or missing weight)'
    return intensity,basis

def treatments(labels):
    treatment=pd.read_csv(PRIVATE/'current_day1_treatments.csv')
    icu=pd.read_csv(W/'icu/data_private/icu_mortality_analysis.csv',
                    usecols=['dataset','icu_stay_id','icu_status'])
    nodes=[];links=[];summaries=[];cells=[];support=[]
    levels=['Low','Medium','High','Max','Intensity unknown']
    for db in ['mimiciv','mimiciii']:
        daily=pd.read_csv(PRIVATE/f'all_icu_d1d7_daily_{db}.csv',low_memory=False)
        daily=daily[daily.day.eq(1)]
        d=labels[labels.dataset.eq(db)][['icu_stay_id','gbmtm_class']].merge(daily,on='icu_stay_id',validate='one_to_one')
        d=d.merge(treatment[treatment.dataset.eq(db)].drop(columns='dataset'),on='icu_stay_id',validate='one_to_one')
        d=d.merge(icu[icu.dataset.eq(db)].drop(columns='dataset'),on='icu_stay_id',validate='one_to_one')
        intensity,fluid_basis=assign_intensity(d)
        valid=intensity.ne('Intensity unknown')
        for c,g in d.groupby('gbmtm_class'):
            for name,values in [('Mechanical ventilation',g.documented_mv_24h.eq(1)),
                ('CRRT',g.documented_crrt_24h.eq(1)),('NEE >0.3 µg/kg/min',g.nee_max_24h.gt(.3)),
                ('Net balance >30 mL/kg',g.fluid_balance_ml_kg_24h.gt(30)),
                ('Absolute fluid cutoff used',fluid_basis.loc[g.index].eq('2,000 mL (invalid or missing weight)'))]:
                support.append(dict(dataset=db,trajectory_class=int(c),component=name,n=int(values.sum()),
                    assigned_n=len(g),uses_outcomes=False))
        for h in [28,'ICU']:
            # Start with the complete assigned cohort. Unobserved intensity or
            # final vital status remains an explicit destination.
            outcome=pd.Series('unknown',index=d.index)
            if h=='ICU':
                outcome[d.icu_status.eq(1)]='icu_death'
                outcome[d.icu_status.eq(2)]='live_exit'
            else:
                outcome[d[f'event_{h}d'].eq(1)]='Died'
                outcome[d[f'event_{h}d'].eq(0)&d[f'followup_{h}d_days'].ge(h)]='Alive'
            paths={0:'Class '+d.gbmtm_class.astype(str),1:intensity,2:outcome}
            nn,ll=aggregate_paths(paths,dict(dataset=db,flow='treatment_intensity',horizon=h))
            nodes.append(nn);links.append(ll)
            summaries.append(dict(dataset=db,exposure='treatment_intensity',horizon=h,full_assigned_n=len(d),
                eligible_n=len(d),excluded_n=0,intensity_known_n=int(valid.sum()),intensity_unknown_n=int((~valid).sum()),
                deaths=int(outcome.isin(['Died','icu_death']).sum()),
                alive=int(outcome.isin(['Alive','live_exit']).sum()),unknown=int(outcome.eq('unknown').sum())))
            cells.extend(pd.DataFrame(dict(trajectory_class=d.gbmtm_class,exposure_level=intensity,outcome=outcome)).groupby(
                ['trajectory_class','exposure_level','outcome']).size().reset_index(name='n').assign(
                    exposure='treatment_intensity',horizon=h,dataset=db).to_dict('records'))
    publish(pd.concat(nodes,ignore_index=True),'supp_treatment_flow_nodes')
    publish(pd.concat(links,ignore_index=True),'supp_treatment_flow_links')
    publish(pd.DataFrame(summaries),'supp_treatment_flow_summary')
    publish(pd.DataFrame(cells),'supp_treatment_flow_cells')
    publish(pd.DataFrame(support),'supp_treatment_intensity_support')
    # Replace the superseded three-bin candidate grid by the actual binary
    # class-by-arm support of the current Figure S5 models.
    publish(pd.read_csv(S/'supp_treatment_support_cells.csv'),'supp_treatment_cutpoint_support')
    distributions=[]
    for db,g in treatment.groupby('dataset'):
        for variable in ['nee_max_24h','fluid_balance_ml_24h','fluid_balance_ml_kg_24h']:
            x=g[variable].dropna();q=x.quantile([.05,.25,.5,.75,.95])
            distributions.append(dict(dataset=db,exposure=variable,n=len(x),
                **{f'p{int(k*100):02d}':v for k,v in q.items()}))
    publish(pd.DataFrame(distributions),'supp_treatment_distributions')
    return dict(treatment_intensity=dict(label='First-day treatment intensity',levels=levels,
        components=['Recorded invasive mechanical ventilation','Recorded CRRT','Maximum NEE','Net fluid balance'],
        hierarchy_order=['Max','High','Medium','Low'],
        hierarchy={'Max':'CRRT or invasive MV with maximum NEE >0.3 µg/kg/min',
            'High':'Below Max: NEE >0.3 or invasive MV with 0<NEE≤0.3 µg/kg/min',
            'Medium':'Below Max/High: invasive MV with NEE=0, isolated 0<NEE≤0.3, or high net fluid balance',
            'Low':'No recorded MV/CRRT, NEE=0, and low net fluid balance'},
        fluid_rule='High if >30 mL/kg with valid 20–300 kg weight; if weight invalid/missing use >2,000 mL; low includes equality',
        unknown='Only no recorded MV/CRRT/NEE with missing fluid documentation, preventing Low/Medium separation',
        selection='User-confirmed organ-support-first descriptive hierarchy; identical in both databases; no outcome optimisation'))

def main():
    labels=pd.read_csv(PRIVATE/'all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    labels=labels[labels.selected_k.eq(3)]
    definitions=treatments(labels)
    disposition(labels)
    record=dict(reference=dict(title='Blood pressure response index trajectories identify distinct hemodynamic phenotypes and predict mortality in septic shock: a two-database retrospective cohort study',
        doi='10.1186/s13054-026-06059-w',url='https://pmc.ncbi.nlm.nih.gov/articles/PMC13154610/',figures=[5,6]),
        fixed_classes=True,classes=[1,2,3],display_days=[1,3,5,7],
        disposition_definition='Fixed class while alive and ICU-resident at each checkpoint; live ICU exit and exact ICU death are absorbing; uncertain date-only exit outcome remains unknown',
        treatment_scope='Both complete assigned cohorts; first-day [0,24 h) treatment intensity; admission Day 28 all-cause death and complete index-ICU outcomes',
        cutpoint_selection='Fixed study-defined clinical categories; no outcome-based optimisation',
        definitions=definitions,zero_semantics='Zero NEE/MV/CRRT denotes no recorded receipt. Intensity is unknown only with no recorded MV/CRRT/NEE and missing fluid documentation; organ support or positive NEE determines intensity despite missing fluid',
        flow_coloring='Persistent source Class 1/2/3; downstream node fills identify exposure/outcome',
        node_percentages='100 × node count / panel cohort_stays, rounded to one decimal; denominator fixed across all stages within each panel',
        models_refitted=False,patient_level_outputs=False)
    target=S/'supplementary_flow_definitions.json';stage=target.with_suffix('.json.tmp')
    stage.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n');stage.replace(target)
    print(pd.read_csv(S/'supp_treatment_flow_summary.csv').to_string(index=False))
    print(pd.read_csv(S/'supp_disposition_summary.csv').to_string(index=False))

if __name__=='__main__':main()
