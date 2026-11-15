"""Check landmark support and generate current risk-set counts."""
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control
W=Path(__file__).resolve().parents[1];A=W/'aggregate'
manifest=pd.read_csv(A/'model_manifest.csv')
assert set(manifest.loc[manifest.primary,'remaining_horizon_days'])=={21,27,83,89}
labels=pd.read_csv(W/'data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
flow=[]
for db in ['mimiciv','mimiciii']:
    d=pd.read_csv(W/'data_private'/f'all_icu_d1d7_daily_{db}.csv',low_memory=False).query('day==1')
    assert d.subject_id.is_unique and d.icu_stay_id.is_unique
    t=d[d.icu_stay_id.isin(labels.loc[labels.dataset.eq(db),'icu_stay_id'])]
    for stage,x in [('cached_adult_first_icu',d),('static_alive_24h',d[d.followup_90d_days>1]),('fixed_k3_assigned',t),('trajectory_alive_icu_resident_day7',t[(t.icu_los_days>=7)&(t.mortality_7d==0)&(t.followup_90d_days>7)])]:
        flow.append(dict(dataset=db,stage=stage,n=len(x),events_28=int(x.event_28d.sum()),events_90=int(x.event_90d.sum())))
pd.DataFrame(flow).to_csv(A/'landmark_population_flow.csv',index=False)
bh=pd.read_csv(A/'all_icu_cox_gbmtm_subgroup_interactions.csv')
for _,g in bh.groupby(['dataset','analysis','horizon']):
    assert len(g)==7
    assert np.allclose(false_discovery_control(g.interaction_p,method='bh'),g.interaction_p_bh)
print('Landmark and BH checks passed')
