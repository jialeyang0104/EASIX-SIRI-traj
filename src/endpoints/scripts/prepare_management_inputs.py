"""Refresh current parent treatment and assigned six-hour inputs (read-only SQL)."""
import pandas as pd
import numpy as np
import psycopg2
from psycopg2.extras import execute_values
import extract_supplementary_treatments as e


def query(db,port,assignments,body,name):
    with psycopg2.connect(dbname=db,port=port) as conn:
        conn.set_session(readonly=True)
        with conn.cursor() as cur:
            rows=execute_values(cur,'WITH assigned(icu_stay_id,trajectory_class) AS (VALUES %s), '+body.replace('%','%%'),list(assignments.itertuples(index=False,name=None)),template='(%s::bigint,%s::integer)',page_size=100000,fetch=True)
            return pd.DataFrame(rows,columns=[x.name for x in cur.description])
e.query=query


def interval_extras(db,port,assigned):
    iv=db=='mimiciv';key='stay_id' if iv else 'icustay_id';schema='mimiciv_derived' if iv else 'mimiciii'
    icu='mimiciv_icu.icustays' if iv else 'mimiciii.icustays'
    mv='ventilation' if iv else 'ventilation_durations';cr='crrt' if iv else 'crrt_durations'
    mvactive="v.ventilation_status='InvasiveVent'" if iv else 'TRUE'
    crtime='c.charttime>=i.t0 AND c.charttime<i.t1' if iv else 'c.starttime<i.t1 AND c.endtime>i.t0'
    cra='(c.system_active=1 OR (c.system_active IS NULL AND c.blood_flow>0))' if iv else 'TRUE'
    body=f'''intervals AS MATERIALIZED (SELECT a.icu_stay_id,k AS interval_index,i.intime,i.intime+(k-1)*INTERVAL '6 hours' t0,i.intime+k*INTERVAL '6 hours' t1 FROM assigned a JOIN {icu} i ON i.{key}=a.icu_stay_id CROSS JOIN generate_series(1,4) k),
    mv AS (SELECT i.icu_stay_id,i.interval_index,COALESCE(MAX(({mvactive})::integer) FILTER(WHERE v.{key} IS NOT NULL),0) AS mv FROM intervals i LEFT JOIN {schema}.{mv} v ON v.{key}=i.icu_stay_id AND v.starttime<i.t1 AND v.endtime>i.t0 GROUP BY 1,2),
    cr AS (SELECT i.icu_stay_id,i.interval_index,COALESCE(MAX(({cra})::integer) FILTER(WHERE c.{key} IS NOT NULL),0) AS crrt FROM intervals i LEFT JOIN {schema}.{cr} c ON c.{key}=i.icu_stay_id AND {crtime} GROUP BY 1,2)
    SELECT i.icu_stay_id,i.interval_index,mv.mv,cr.crrt'''
    if iv:
        body+=''',f.input_events,f.input_ml,o.output_events,o.output_ml FROM intervals i
        LEFT JOIN (SELECT i.icu_stay_id,i.interval_index,COUNT(f.amount) input_events,SUM(CASE WHEN f.endtime>f.starttime THEN f.amount*EXTRACT(EPOCH FROM(LEAST(f.endtime,i.t1)-GREATEST(f.starttime,i.t0)))/EXTRACT(EPOCH FROM(f.endtime-f.starttime)) ELSE f.amount END) input_ml FROM intervals i LEFT JOIN mimiciv_icu.inputevents f ON f.stay_id=i.icu_stay_id AND f.starttime<i.t1 AND ((f.endtime>f.starttime AND f.endtime>i.t0) OR (f.endtime=f.starttime AND f.starttime>=i.t0)) AND LOWER(f.amountuom)='ml' AND f.amount>0 AND COALESCE(f.statusdescription,'') NOT IN ('Rewritten','Cancelled') GROUP BY 1,2) f USING(icu_stay_id,interval_index)
        LEFT JOIN (SELECT i.icu_stay_id,i.interval_index,COUNT(o.value) output_events,SUM(o.value) output_ml FROM intervals i LEFT JOIN mimiciv_icu.outputevents o ON o.stay_id=i.icu_stay_id AND o.charttime>=i.t0 AND o.charttime<i.t1 AND o.value>0 AND LOWER(o.valueuom)='ml' GROUP BY 1,2) o USING(icu_stay_id,interval_index)'''
    else:body+=' FROM intervals i'
    body+=' LEFT JOIN mv USING(icu_stay_id,interval_index) LEFT JOIN cr USING(icu_stay_id,interval_index)'
    return query(db,port,assigned,body,'unused')


def main():
    labels=pd.read_csv(e.PRIVATE/'all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    frames=[]
    for db,port in [('mimiciv',5432),('mimiciii',5433)]:
        daily=pd.read_csv(e.PRIVATE/f'all_icu_d1d7_daily_{db}.csv',usecols=['icu_stay_id','day'])
        assigned=daily[daily.day.eq(1)][['icu_stay_id']].copy();assigned['trajectory_class']=0
        print(db,'parent',len(assigned),flush=True)
        s=e.support(db,port,assigned)
        if db=='mimiciv':
            t=pd.read_csv(e.P/'data/private/all_icu_day1_treatments_mimiciv.csv')
        else:
            prof=e.mimiciii_profiles(assigned)
            g=prof.groupby('icu_stay_id');t=g.agg(weight_raw_kg=('weight_raw_kg','first'),weight_kg=('weight_kg','first'),nee_max_24h=('nee_max','max'),nee_any_24h=('nee_any','max'),nee_unquantified=('nee_unquantified','max'),n_input_events=('input_events','sum'),n_output_events=('output_events','sum'),fluid_input_ml_24h=('input_ml',lambda x:x.sum(min_count=1)),fluid_output_ml_24h=('output_ml',lambda x:x.sum(min_count=1))).reset_index()
            t.loc[t.nee_unquantified.eq(1),'nee_max_24h']=np.nan
            t['input_documented']=t.n_input_events.gt(0).astype(int);t['output_documented']=t.n_output_events.gt(0).astype(int)
            t['fluid_balance_ml_24h']=t.fluid_input_ml_24h-t.fluid_output_ml_24h;t['fluid_balance_ml_kg_24h']=t.fluid_balance_ml_24h/t.weight_kg
        t['dataset']=db;t=assigned[['icu_stay_id']].merge(t,on='icu_stay_id',validate='one_to_one').merge(s,on='icu_stay_id',validate='one_to_one');frames.append(t)
        lab=labels[labels.dataset.eq(db)][['icu_stay_id','gbmtm_class']].astype(int)
        extras=interval_extras(db,port,lab)
        if db=='mimiciii':
            p=prof[prof.icu_stay_id.isin(lab.icu_stay_id)].drop(columns='trajectory_class').merge(lab.rename(columns={'gbmtm_class':'trajectory_class'}),on='icu_stay_id')
        else:
            p=pd.read_csv(e.P/'data/private/mimiciv_same_interval_profiles_patient_level.csv');p=p[p.interval_index.le(4)&p.icu_stay_id.isin(lab.icu_stay_id)]
        p=p.merge(extras,on=['icu_stay_id','interval_index'],validate='one_to_one');e.publish(p,f'{db}_current_six_hour_management.csv')
        print(db,'done',len(t),len(p),flush=True)
    e.publish(pd.concat(frames,ignore_index=True),'current_parent_day1_treatments.csv')
if __name__=='__main__':main()
