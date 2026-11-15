"""Current-cohort first-day treatment sources for supplementary Figures S4–S7.

Read-only database queries use the current fixed assignments, without the old
MIMIC-III admission-level chart-event restriction. Patient rows stay private.
"""
from pathlib import Path
from datetime import datetime, timezone
import json

import numpy as np
import pandas as pd
import psycopg2
from psycopg2.extras import execute_values

W = Path(__file__).resolve().parents[1]
P = W.parents[1]
OLD = P
PRIVATE = W/'data_private'


def query(db, port, assignments, body, name):
    sql = '''WITH assigned(icu_stay_id,trajectory_class) AS (VALUES %s), ''' + body.replace('%','%%')
    (W/'sql'/name).write_text(sql + '\n')
    with psycopg2.connect(dbname=db, port=port) as connection:
        connection.set_session(readonly=True)
        with connection.cursor() as cursor:
            rows = execute_values(cursor, sql, list(assignments.itertuples(index=False, name=None)),
                                  template='(%s::bigint,%s::integer)', page_size=5000, fetch=True)
            return pd.DataFrame(rows, columns=[c.name for c in cursor.description])


def support(db, port, assignments):
    iv = db == 'mimiciv'
    key = 'stay_id' if iv else 'icustay_id'
    icu = 'mimiciv_icu.icustays' if iv else 'mimiciii.icustays'
    schema = 'mimiciv_derived' if iv else 'mimiciii'
    ventilation = 'ventilation' if iv else 'ventilation_durations'
    mv = "v.ventilation_status='InvasiveVent'" if iv else 'TRUE'
    cr_table = 'crrt' if iv else 'crrt_durations'
    cr_window = ("c.charttime>=e.intime AND c.charttime<e.intime+INTERVAL '24 hours'"
                 if iv else "c.starttime<e.intime+INTERVAL '24 hours' AND c.endtime>e.intime")
    cr_active = '(c.system_active=1 OR (c.system_active IS NULL AND c.blood_flow>0))' if iv else 'TRUE'
    body = f'''cohort AS MATERIALIZED (SELECT i.*,a.trajectory_class FROM assigned a JOIN {icu} i
        ON i.{key}=a.icu_stay_id), mv AS (
            SELECT e.{key},COUNT(v.{key}) AS source_rows,
                COALESCE(MAX(({mv})::integer) FILTER(WHERE v.{key} IS NOT NULL),0) AS active
            FROM cohort e LEFT JOIN {schema}.{ventilation} v ON v.{key}=e.{key}
                AND v.starttime<e.intime+INTERVAL '24 hours' AND v.endtime>e.intime GROUP BY e.{key}
        ), cr AS (
            SELECT e.{key},COUNT(c.{key}) AS source_rows,
                COALESCE(MAX(({cr_active})::integer) FILTER(WHERE c.{key} IS NOT NULL),0) AS active
            FROM cohort e LEFT JOIN {schema}.{cr_table} c ON c.{key}=e.{key} AND {cr_window} GROUP BY e.{key}
        )
        SELECT e.{key} AS icu_stay_id, o.mechvent_score AS oasis_mv_points,
            o.oasis-o.mechvent_score AS oasis_without_mv,
            mv.source_rows AS mv_source_rows,mv.active AS documented_mv_24h,
            cr.source_rows AS crrt_source_rows,cr.active AS documented_crrt_24h
        FROM cohort e LEFT JOIN {schema}.oasis o ON o.{key}=e.{key}
        LEFT JOIN mv ON mv.{key}=e.{key} LEFT JOIN cr ON cr.{key}=e.{key} ORDER BY e.{key}'''
    return query(db, port, assignments, body, f'{db}_day1_documented_support.sql')


def mimiciii_profiles(assignments):
    # Reuse the locally audited CareVue unit corrections; the cohort, fluid
    # mapping and chart-event eligibility are replaced here, not in the archive.
    old_sql = (OLD/'sql/07_mimiciii_all_icu_day1_treatments.sql').read_text()
    weight_and_units = '''weight AS MATERIALIZED (
        SELECT DISTINCT ON(e.icustay_id) e.icustay_id,w.weight AS weight_raw_kg,
            CASE WHEN w.weight BETWEEN 20 AND 300 THEN w.weight END AS weight_kg
        FROM eligible e LEFT JOIN mimiciii.weight_durations w ON w.icustay_id=e.icustay_id
            AND w.starttime<e.intime+INTERVAL '24 hours' AND w.endtime>e.intime
        ORDER BY e.icustay_id,CASE WHEN w.starttime<=e.intime AND w.endtime>e.intime THEN 0 ELSE 1 END,w.starttime
    ), ''' + old_sql[old_sql.index('phenylephrine_cv_mcg_min AS MATERIALIZED'):old_sql.index('input_union AS MATERIALIZED')]
    component = old_sql[old_sql.index('component_events AS MATERIALIZED'):old_sql.index('nee_boundaries AS (')]
    body = '''eligible AS MATERIALIZED (
        SELECT i.*,a.trajectory_class FROM assigned a JOIN mimiciii.icustays i
        ON i.icustay_id=a.icu_stay_id
    ), ''' + weight_and_units + component + '''
    bounds AS (
        SELECT icustay_id,intime AS t FROM eligible UNION
        SELECT icustay_id,intime+INTERVAL '24 hours' FROM eligible UNION
        SELECT n.icustay_id,GREATEST(n.starttime,e.intime) FROM component_events n JOIN eligible e USING(icustay_id)
        WHERE n.starttime<e.intime+INTERVAL '24 hours' AND n.endtime>e.intime UNION
        SELECT n.icustay_id,LEAST(n.endtime,e.intime+INTERVAL '24 hours') FROM component_events n JOIN eligible e USING(icustay_id)
        WHERE n.starttime<e.intime+INTERVAL '24 hours' AND n.endtime>e.intime
    ), segments AS (
        SELECT icustay_id,t AS starttime,LEAD(t) OVER(PARTITION BY icustay_id ORDER BY t) AS endtime FROM bounds
    ), per_drug AS (
        SELECT s.icustay_id,s.starttime,s.endtime,n.drug,MAX(n.rate*n.factor) AS dose,
            BOOL_OR(n.drug IS NOT NULL AND n.rate IS NULL) AS unquantified
        FROM segments s LEFT JOIN component_events n ON n.icustay_id=s.icustay_id
        AND n.starttime<=s.starttime AND n.endtime>s.starttime
        WHERE s.endtime>s.starttime GROUP BY s.icustay_id,s.starttime,s.endtime,n.drug
    ), nee AS (
        SELECT icustay_id,starttime,endtime,SUM(dose) AS dose,
            BOOL_OR(unquantified) AS unquantified FROM per_drug GROUP BY 1,2,3
    ), intervals AS (
        SELECT e.*,k AS interval_index,e.intime+(k-1)*INTERVAL '6 hours' AS t0,
            e.intime+k*INTERVAL '6 hours' AS t1,w.weight_raw_kg,w.weight_kg
        FROM eligible e CROSS JOIN generate_series(1,4) k
        LEFT JOIN weight w ON w.icustay_id=e.icustay_id
    ), input_events AS MATERIALIZED (
        SELECT im.icustay_id,im.starttime,COALESCE(im.endtime,im.starttime) AS endtime,im.amount,
            (LOWER(COALESCE(im.ordercategoryname,'')) SIMILAR TO '(02-%|03-%|04-%)'
                AND (LOWER(COALESCE(di.label,'')) ~ 'albumin|hetastarch|dextran|saline|nacl|ringer|(^|[^a-z])lr([^a-z]|$)|lactated|dextrose|d5|sterile water|crystalloid'
                     OR LOWER(COALESCE(di.category,''))='fluids/intake')
                AND LOWER(COALESCE(di.label,'')) NOT IN ('piggyback','solution','multivitamins')
                AND LOWER(COALESCE(di.label,'')) !~ 'tube feed|replete|nepro|osmolite|jevity|promote|glucerna|nutren|packed red|platelet|fresh frozen|cryoprecip') AS audited,
            'MetaVision'::text AS source
        FROM eligible e JOIN mimiciii.inputevents_mv im USING(icustay_id)
        JOIN mimiciii.d_items di USING(itemid)
        WHERE im.starttime<e.intime+INTERVAL '24 hours' AND COALESCE(im.endtime,im.starttime)>=e.intime
            AND LOWER(im.amountuom)='ml' AND im.amount>0 AND im.cancelreason=0
            AND COALESCE(im.statusdescription,'') NOT IN ('Rewritten','Cancelled')
        UNION ALL
        SELECT cv.icustay_id,cv.charttime,cv.charttime,cv.amount,
            cv.itemid IN (30008,30009,30018,30020,30021,30143,30159,30160,30161,30176,30181,30352,30353,
                40548,42832,43237,43353,44203,44440,44952,45403,45532,46564),'CareVue'
        FROM eligible e JOIN mimiciii.inputevents_cv cv USING(icustay_id)
        WHERE cv.charttime>=e.intime AND cv.charttime<e.intime+INTERVAL '24 hours'
            AND LOWER(cv.amountuom)='ml' AND cv.amount>0
    ), input_interval AS (
        SELECT i.icustay_id,i.interval_index,COUNT(f.amount) AS input_events,
            SUM(CASE WHEN f.endtime>f.starttime THEN f.amount*EXTRACT(EPOCH FROM(LEAST(f.endtime,i.t1)-GREATEST(f.starttime,i.t0)))
                /EXTRACT(EPOCH FROM(f.endtime-f.starttime)) ELSE f.amount END) AS input_ml,
            COALESCE(SUM(CASE WHEN f.endtime>f.starttime THEN f.amount*EXTRACT(EPOCH FROM(LEAST(f.endtime,i.t1)-GREATEST(f.starttime,i.t0)))
                /EXTRACT(EPOCH FROM(f.endtime-f.starttime)) ELSE f.amount END) FILTER(WHERE f.audited),0) AS audited_iv_ml,
            COUNT(f.amount) FILTER(WHERE f.source='MetaVision') AS mv_input_events,
            COUNT(f.amount) FILTER(WHERE f.source='CareVue') AS cv_input_events
        FROM intervals i LEFT JOIN input_events f ON f.icustay_id=i.icustay_id AND f.starttime<i.t1
            AND ((f.endtime>f.starttime AND f.endtime>i.t0) OR (f.endtime=f.starttime AND f.starttime>=i.t0)) GROUP BY 1,2
    ), output_interval AS (
        SELECT i.icustay_id,i.interval_index,COUNT(o.value) AS output_events,SUM(o.value) AS output_ml
        FROM intervals i LEFT JOIN mimiciii.outputevents o ON o.icustay_id=i.icustay_id
            AND o.charttime>=i.t0 AND o.charttime<i.t1 AND LOWER(o.valueuom)='ml' AND o.value>0 AND COALESCE(o.iserror,0)=0 GROUP BY 1,2
    ), nee_interval AS (
        SELECT i.icustay_id,i.interval_index,CASE WHEN BOOL_OR(n.unquantified) THEN NULL ELSE COALESCE(MAX(n.dose),0) END AS nee_max,
            (COALESCE(MAX(n.dose),0)>0 OR COALESCE(BOOL_OR(n.unquantified),FALSE))::integer AS nee_any,
            COALESCE(BOOL_OR(n.unquantified),FALSE)::integer AS nee_unquantified
        FROM intervals i LEFT JOIN nee n ON n.icustay_id=i.icustay_id AND n.starttime<i.t1 AND n.endtime>i.t0 GROUP BY 1,2
    ), vitals AS (
        SELECT i.icustay_id,i.interval_index,COUNT(v.meanbp) FILTER(WHERE v.meanbp BETWEEN 20 AND 200) AS map_observed_n,
            percentile_cont(.5) WITHIN GROUP(ORDER BY v.meanbp) FILTER(WHERE v.meanbp BETWEEN 20 AND 200) AS map_median
        FROM intervals i LEFT JOIN mimiciii.pivoted_vital v ON v.icustay_id=i.icustay_id AND v.charttime>=i.t0 AND v.charttime<i.t1 GROUP BY 1,2
    ), lactate AS (
        SELECT i.icustay_id,i.interval_index,COUNT(l.valuenum) AS lactate_observed_n,
            percentile_cont(.5) WITHIN GROUP(ORDER BY l.valuenum) AS lactate_median
        FROM intervals i LEFT JOIN mimiciii.labevents l ON l.hadm_id=i.hadm_id AND l.charttime>=i.t0 AND l.charttime<i.t1
            AND l.itemid=50813 AND l.valuenum>0 GROUP BY 1,2
    ), urine AS (
        SELECT i.icustay_id,i.interval_index,COUNT(u.urineoutput) AS urine_observed_n,SUM(u.urineoutput) AS urine_output_ml
        FROM intervals i LEFT JOIN mimiciii.pivoted_uo u ON u.icustay_id=i.icustay_id AND u.charttime>=i.t0 AND u.charttime<i.t1
            AND u.urineoutput>=0 GROUP BY 1,2
    )
    SELECT i.icustay_id AS icu_stay_id,i.trajectory_class,i.interval_index,i.weight_raw_kg,i.weight_kg,
        1 AS at_risk_start,6 AS observed_window_hours,n.nee_max,n.nee_any,n.nee_unquantified,
        f.input_events,f.input_ml,f.audited_iv_ml,
        CASE WHEN f.audited_iv_ml=0 THEN 0 ELSE f.audited_iv_ml/i.weight_kg END AS audited_iv_ml_kg,
        f.mv_input_events,f.cv_input_events,o.output_events,o.output_ml,
        v.map_observed_n,v.map_median,l.lactate_observed_n,l.lactate_median,u.urine_observed_n,
        u.urine_output_ml,u.urine_output_ml/i.weight_kg/6.0 AS urine_output_ml_kg_h
    FROM intervals i LEFT JOIN nee_interval n USING(icustay_id,interval_index)
    LEFT JOIN input_interval f USING(icustay_id,interval_index) LEFT JOIN output_interval o USING(icustay_id,interval_index)
    LEFT JOIN vitals v USING(icustay_id,interval_index) LEFT JOIN lactate l USING(icustay_id,interval_index)
    LEFT JOIN urine u USING(icustay_id,interval_index) ORDER BY i.icustay_id,i.interval_index'''
    return query('mimiciii',5433,assignments,body,'mimiciii_six_hour_treatments.sql')


def publish(frame, name):
    target=PRIVATE/name; stage=target.with_suffix('.csv.tmp')
    frame.to_csv(stage,index=False);stage.replace(target)


def main():
    labels=pd.read_csv(PRIVATE/'all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
    labels=labels[labels.selected_k.eq(3)]
    frames=[];audit=[]
    for db,port in [('mimiciv',5432),('mimiciii',5433)]:
        assigned=labels[labels.dataset.eq(db)][['icu_stay_id','gbmtm_class']].astype(int)
        s=support(db,port,assigned)
        if db=='mimiciv':
            t=pd.read_csv(OLD/'data/private/all_icu_day1_treatments_mimiciv.csv')
            t=assigned[['icu_stay_id']].merge(t,on='icu_stay_id',validate='one_to_one')
        else:
            prof=mimiciii_profiles(assigned)
            assert len(prof)==4*len(assigned) and not prof.duplicated(['icu_stay_id','interval_index']).any()
            publish(prof,'mimiciii_six_hour_management.csv')
            groups=prof.groupby('icu_stay_id')
            t=groups.agg(weight_raw_kg=('weight_raw_kg','first'),weight_kg=('weight_kg','first'),
                nee_max_24h=('nee_max','max'),nee_any_24h=('nee_any','max'),nee_unquantified=('nee_unquantified','max'),
                n_input_events=('input_events','sum'),n_output_events=('output_events','sum'),
                fluid_input_ml_24h=('input_ml',lambda v:v.sum(min_count=1)),
                fluid_output_ml_24h=('output_ml',lambda v:v.sum(min_count=1)),
                mv_input_events=('mv_input_events','sum'),cv_input_events=('cv_input_events','sum')).reset_index()
            t.loc[t.nee_unquantified.eq(1),'nee_max_24h']=np.nan
            t['input_documented']=t.n_input_events.gt(0).astype(int)
            t['output_documented']=t.n_output_events.gt(0).astype(int)
            t['fluid_balance_ml_24h']=t.fluid_input_ml_24h-t.fluid_output_ml_24h
            t['fluid_balance_ml_kg_24h']=t.fluid_balance_ml_24h/t.weight_kg
            t['dataset']=db
        columns=['dataset','icu_stay_id','weight_raw_kg','weight_kg','n_input_events','n_output_events',
            'input_documented','output_documented','fluid_input_ml_24h','fluid_output_ml_24h',
            'fluid_balance_ml_24h','fluid_balance_ml_kg_24h','nee_any_24h','nee_max_24h']
        t=t[columns].merge(s,on='icu_stay_id',validate='one_to_one')
        assert len(t)==len(assigned) and t.icu_stay_id.is_unique
        frames.append(t)
        audit.append(dict(dataset=db,assigned_n=len(assigned),extracted_n=len(t),
            input_documented=int(t.input_documented.sum()),output_documented=int(t.output_documented.sum()),
            io_documented=int((t.input_documented.eq(1)&t.output_documented.eq(1)).sum()),
            weight_valid=int(t.weight_kg.between(20,300).sum()),
            recorded_mv=int(t.documented_mv_24h.sum()),recorded_crrt=int(t.documented_crrt_24h.sum()),
            mv_source_no_row=int(t.mv_source_rows.eq(0).sum()),crrt_source_no_row=int(t.crrt_source_rows.eq(0).sum()),
            nee_missing=int(t.nee_max_24h.isna().sum())))
        print(db,audit[-1],flush=True)
    publish(pd.concat(frames,ignore_index=True),'current_day1_treatments.csv')
    pd.DataFrame(audit).to_csv(W/'aggregate/treatment_source_coverage.csv',index=False)
    record=dict(extracted_at=datetime.now(timezone.utc).isoformat(),current_cohort=True,chart_event_flag_required=False,
        iv_treatment_source=str((OLD/'data/private/all_icu_day1_treatments_mimiciv.csv').relative_to(P)),
        iii_treatment_source='Read-only current-cohort query; all positive-mL MetaVision/CareVue input and error-free output',
        cv_amount_rule='Positive recorded amounts retained, including stopped records; stopped indicates infusion cessation, not cancellation',
        nee_formula='NE + E + phenylephrine/10 + dopamine/100 + 2.5*vasopressin(U/min); IV treatment cache additionally includes 10*angiotensinII',
        mechanical_ventilation='Recorded invasive-ventilation interval overlapping [0,24 h); no positive record is no recorded receipt',
        crrt='IV: system_active=1 or missing system_active with blood_flow>0; III: derived CRRT duration overlaps [0,24 h)',
        support_zero_semantics='0=no recorded receipt, not confirmed absence; raw source-row counts retained',
        adjustment='OASIS minus its mechanical-ventilation points, consistently across all S5 exposure models',
        fluid_zero_semantics='Missing input or output remains missing, not zero',patient_rows_packaged=False)
    (W/'aggregate/treatment_source_manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')


if __name__=='__main__':
    main()
