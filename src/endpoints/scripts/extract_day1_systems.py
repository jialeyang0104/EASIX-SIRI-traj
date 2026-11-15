"""Strict post-ICU 0–24 h measurement means, harmonised across database versions."""
from pathlib import Path
import subprocess, io, json
import pandas as pd

P=Path(__file__).resolve().parents[3]
W=P/'src/endpoints'
# Same laboratory item identifiers and units verified in both d_labitems tables.
LABS=[
 (50861,'alt_mean','ALT','U/L','Hepatic'),(50878,'ast_mean','AST','U/L','Hepatic'),
 (50863,'alp_mean','ALP','U/L','Hepatic'),(50885,'bilirubin_mean','Total bilirubin','mg/dL','Hepatic'),
 (50862,'albumin_mean','Albumin','g/dL','Hepatic'),(50912,'creatinine_mean','Creatinine','mg/dL','Renal'),
 (51006,'bun_mean','BUN','mg/dL','Renal'),(50954,'ldh_mean','LDH','U/L','Index components'),
 (51301,'wbc_mean','WBC','10^9/L','Haematology'),(51265,'platelet_mean','Platelets','10^9/L','Haematology'),
 (51222,'hemoglobin_mean','Haemoglobin','g/dL','Haematology'),(51221,'hematocrit_mean','Haematocrit','%','Haematology'),
 (51279,'rbc_mean','RBC','10^12/L','Haematology'),(51250,'mcv_mean','MCV','fL','Haematology'),
 (51237,'inr_mean','INR','ratio','Coagulation'),(51274,'pt_mean','PT','s','Coagulation'),
 (51275,'ptt_mean','aPTT','s','Coagulation'),(50983,'sodium_mean','Sodium','mmol/L','Electrolytes'),
 (50971,'potassium_mean','Potassium','mmol/L','Electrolytes'),(50902,'chloride_mean','Chloride','mmol/L','Electrolytes'),
 (50882,'bicarbonate_mean','Bicarbonate','mmol/L','Electrolytes'),(50868,'aniongap_mean','Anion gap','mmol/L','Metabolic'),
 (50893,'calcium_mean','Total calcium','mg/dL','Electrolytes'),(50960,'magnesium_mean','Magnesium','mg/dL','Electrolytes'),
 (50970,'phosphate_mean','Phosphate','mg/dL','Electrolytes'),(50931,'glucose_mean','Glucose','mg/dL','Metabolic'),
 (50813,'lactate_mean','Lactate','mmol/L','Metabolic')]
VITALS=[('heart_rate','heartrate','heart_rate_mean','Heart rate','beats/min','Cardiovascular'),
 ('sbp','sysbp','sbp_mean','SBP','mmHg','Cardiovascular'),('dbp','diasbp','dbp_mean','DBP','mmHg','Cardiovascular'),
 ('mbp','meanbp','map_mean','MAP','mmHg','Cardiovascular'),('resp_rate','resprate','resp_rate_mean','Respiratory rate','breaths/min','Respiratory'),
 ('temperature','tempc','temperature_mean','Temperature','degC','Metabolic'),('spo2','spo2','spo2_mean','SpO2','%','Respiratory')]
BGS=[('ph','ph_mean','Arterial pH','pH','Respiratory'),('po2','pao2_mean','PaO2','mmHg','Respiratory'),
 ('pco2','paco2_mean','PaCO2','mmHg','Respiratory'),('pao2fio2ratio','pf_ratio_mean','PaO2/FiO2','mmHg','Respiratory'),
 ('baseexcess','baseexcess_mean','Base excess','mmol/L','Metabolic')]

def read(db,port,sql):
    r=subprocess.run(['psql','-Xq','-w','-v','ON_ERROR_STOP=1','-p',str(port),'-d',db],input='BEGIN READ ONLY;\nCOPY ('+sql+') TO STDOUT WITH CSV HEADER;\nCOMMIT;',text=True,capture_output=True,check=True)
    return pd.read_csv(io.StringIO(r.stdout))

metadata=[]
for db,port,num,key in [('mimiciv',5432,'01','stay_id'),('mimiciii',5433,'02','icustay_id')]:
    src=(W/'sql'/f'{num}_{db}_all_icu_d1d7_daily.sql').read_text()
    end=src.index('\n),',src.index('eligible AS'))+2
    prefix=src[:end]
    schema='mimiciv_hosp' if db=='mimiciv' else 'mimiciii'
    vital='mimiciv_derived.vitalsign' if db=='mimiciv' else 'mimiciii.pivoted_vital'
    gcs='mimiciv_derived.gcs' if db=='mimiciv' else 'mimiciii.pivoted_gcs'
    bg='mimiciv_derived.bg' if db=='mimiciv' else 'mimiciii.pivoted_bg_art'
    uo='mimiciv_derived.urine_output' if db=='mimiciv' else 'mimiciii.pivoted_uo'
    lcols=',\n'.join(f'AVG(l.valuenum) FILTER (WHERE l.itemid={i} AND l.valuenum {">=0" if var in ["aniongap_mean","bilirubin_mean"] else ">0"}) AS {var}' for i,var,*_ in LABS)
    vcols=',\n'.join(f'AVG(v.{r[0] if db=="mimiciv" else r[1]}) AS {r[2]}' for r in VITALS)
    bcols=',\n'.join(f'AVG(b.{c}) AS {var}' for c,var,*_ in BGS)
    art="AND b.specimen='ART.'" if db=='mimiciv' else ''
    sql=prefix+f''', labs AS (
        SELECT e.{key}, {lcols} FROM eligible e JOIN {schema}.labevents l
        ON l.hadm_id=e.hadm_id AND l.charttime>=e.intime AND l.charttime<e.intime+INTERVAL '24 hours'
        WHERE l.itemid IN ({','.join(str(r[0]) for r in LABS)}) GROUP BY e.{key}
    ), vitals AS (
        SELECT e.{key}, {vcols} FROM eligible e JOIN {vital} v
        ON v.{key}=e.{key} AND v.charttime>=e.intime AND v.charttime<e.intime+INTERVAL '24 hours' GROUP BY e.{key}
    ), neurologic AS (
        SELECT e.{key}, AVG(g.gcs) AS gcs_mean FROM eligible e JOIN {gcs} g
        ON g.{key}=e.{key} AND g.charttime>=e.intime AND g.charttime<e.intime+INTERVAL '24 hours' GROUP BY e.{key}
    ), bloodgas AS (
        SELECT e.{key}, {bcols} FROM eligible e JOIN {bg} b
        ON b.{'hadm_id' if db=='mimiciv' else key}=e.{'hadm_id' if db=='mimiciv' else key}
        AND b.charttime>=e.intime AND b.charttime<e.intime+INTERVAL '24 hours' {art} GROUP BY e.{key}
    ), urine AS (
        SELECT e.{key}, SUM(u.urineoutput)/24.0 AS urine_hourly_mean FROM eligible e JOIN {uo} u
        ON u.{key}=e.{key} AND u.charttime>=e.intime AND u.charttime<e.intime+INTERVAL '24 hours'
        WHERE u.urineoutput>=0 GROUP BY e.{key}
    ) SELECT e.{key} AS icu_stay_id,
        {','.join('l.'+r[1] for r in LABS)}, {','.join('v.'+r[2] for r in VITALS)},
        g.gcs_mean, {','.join('b.'+r[1] for r in BGS)}, u.urine_hourly_mean
        FROM eligible e LEFT JOIN labs l USING ({key}) LEFT JOIN vitals v USING ({key})
        LEFT JOIN neurologic g USING ({key}) LEFT JOIN bloodgas b USING ({key}) LEFT JOIN urine u USING ({key}) ORDER BY e.{key}'''
    (W/'sql'/f'{db}_day1_system_means.sql').write_text('BEGIN READ ONLY;\nCOPY ('+sql+') TO STDOUT WITH CSV HEADER;\nCOMMIT;\n')
    d=read(db,port,sql)
    assert d.icu_stay_id.is_unique
    assert len(d)==len(pd.read_csv(W/'data_private'/f'{db}_raw_endpoints.csv'))
    path=W/'data_private'/f'{db}_day1_system_means.csv'
    d.to_csv(path.with_suffix('.tmp'),index=False);path.with_suffix('.tmp').replace(path)
    for i,var,label,unit,system in LABS:
        metadata.append(dict(dataset=db,variable=var,label=label,unit=unit,system=system,
            source=f'{schema}.labevents',source_itemid=str(i),window='[ICU admission, ICU admission +24 h)',aggregation='arithmetic mean of available positive measurements; zero allowed for bilirubin and anion gap',index_component=var in ['ldh_mean','creatinine_mean','platelet_mean','wbc_mean']))
    for c1,c2,var,label,unit,system in VITALS:
        metadata.append(dict(dataset=db,variable=var,label=label,unit=unit,system=system,source=vital,source_itemid=c1 if db=='mimiciv' else c2,window='[0,24) ICU hours',aggregation='arithmetic mean of available concept-cleaned measurements',index_component=False))
    for c,var,label,unit,system in BGS:
        metadata.append(dict(dataset=db,variable=var,label=label,unit=unit,system=system,source=bg,source_itemid=c,window='[0,24) ICU hours',aggregation='arithmetic mean of arterial blood gas measurements',index_component=False))
    metadata.extend([dict(dataset=db,variable='gcs_mean',label='GCS',unit='points',system='Neurologic',source=gcs,source_itemid='gcs',window='[0,24) ICU hours',aggregation='arithmetic mean of concept-derived GCS',index_component=False),dict(dataset=db,variable='urine_hourly_mean',label='Urine output',unit='mL/h',system='Renal',source=uo,source_itemid='urineoutput',window='[0,24) ICU hours',aggregation='sum of documented non-negative outputs / 24 h; no records remains missing',index_component=False)])
    print(db, len(d), 'patients;', d.shape[1]-1, 'mean variables',flush=True)
pd.DataFrame(metadata).to_csv(W/'aggregate/day1_system_variable_dictionary.csv',index=False)
