"""Extract licensed source data in read-only transactions. Run at repository root."""
from pathlib import Path
import io
import subprocess
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
END=ROOT/'src/endpoints'
for path in [ROOT/'data/private',ROOT/'results',ROOT/'results/trajectory/tables',ROOT/'results/trajectory/models',END/'data_private',END/'aggregate',END/'logs',ROOT/'src/classifier/aggregate',ROOT/'src/classifier/data_private',ROOT/'src/icu/aggregate',ROOT/'src/icu/data_private',ROOT/'src/icu/logs']:
    path.mkdir(parents=True,exist_ok=True)
def query(db,port,sql):
    return subprocess.run(['psql','-Xq','-w','-v','ON_ERROR_STOP=1','-p',str(port),'-d',db],input='BEGIN READ ONLY;\nCOPY ('+sql.strip().rstrip(';')+') TO STDOUT WITH CSV HEADER;\nCOMMIT;\n',text=True,capture_output=True,check=True).stdout
for db,port,num,key in [('mimiciv',5432,'01','stay_id'),('mimiciii',5433,'02','icustay_id')]:
    sql=(END/'sql'/f'{num}_{db}_all_icu_d1d7_daily.sql').read_text()
    daily=query(db,port,sql)
    (END/'data_private'/f'all_icu_d1d7_daily_{db}.csv').write_text(daily)
    (ROOT/'results/trajectory/tables'/f'all_icu_d1d7_daily_{db}.csv').write_text(daily)
    boundary=sql.index('\n),',sql.index('eligible AS'))+2
    endpoint=sql[:boundary]+f" SELECT {key} AS icu_stay_id,subject_id,hadm_id AS encounter_id,intime,outtime,admittime,dischtime,deathtime,dod,analysis_deathtime AS corrected_death FROM eligible"
    (END/'data_private'/f'{db}_raw_endpoints.csv').write_text(query(db,port,endpoint))
    ranked=sql[:sql.index('\n),',sql.index('WITH icu_ranked'))+2]
    stages=[('All ICU stays','TRUE'),('First ICU stay per patient','patient_icu_rank=1'),('Age >=18 years','patient_icu_rank=1 AND age>=18'),('ICU length of stay >=24 h',"patient_icu_rank=1 AND age>=18 AND outtime-intime>=INTERVAL '24 hours'"),('Hospital length of stay >=24 h',"patient_icu_rank=1 AND age>=18 AND outtime-intime>=INTERVAL '24 hours' AND dischtime-admittime>=INTERVAL '24 hours'")]
    q=ranked+' SELECT '+','.join(f'COUNT(*) FILTER (WHERE {c}) AS n{i}' for i,(_,c) in enumerate(stages))+' FROM icu_ranked'
    counts=pd.read_csv(io.StringIO(query(db,port,q))).iloc[0].astype(int).tolist()
    pd.DataFrame([dict(dataset=db,step_order=i+1,step=label,n=counts[i],excluded_from_previous=None if i==0 else counts[i-1]-counts[i]) for i,(label,_) in enumerate(stages)]).to_csv(END/'aggregate'/f'all_icu_flow_counts_{db}.csv',index=False)
sql=(ROOT/'sql/06_mimiciv_all_icu_day1_treatments.sql').read_text()
(ROOT/'data/private/all_icu_day1_treatments_mimiciv.csv').write_text(query('mimiciv',5432,sql))
