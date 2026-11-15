"""Train-only three-way consensus selection; untouched IV/III validation and SHAP."""
from pathlib import Path
import json, sys, joblib
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.base import clone
from sklearn.preprocessing import QuantileTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.metrics import (roc_auc_score,average_precision_score,accuracy_score,balanced_accuracy_score,
    f1_score,log_loss,roc_curve,confusion_matrix)
from xgboost import XGBClassifier
import xgboost, sklearn, shap

P=Path(__file__).resolve().parents[2]
W=P/'src'
H=W/'classifier';A=H/'aggregate';PRIVATE=H/'data_private'
A.mkdir(exist_ok=True);PRIVATE.mkdir(exist_ok=True)
SEED=20261111
rng=np.random.default_rng(SEED)
labels=pd.read_csv(W/'endpoints/data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv')
meta=pd.read_csv(W/'endpoints/aggregate/day1_system_variable_dictionary.csv').drop_duplicates('variable')
base=['age','male','mechanical_ventilation_day1','vasopressor_day1','crrt_day1','easix_day','siri_day','ldh_u_l','creatinine_mg_dl','platelets_109_l','wbc_109_l','neutrophils_harmonised','monocytes_harmonised','lymphocytes_harmonised']
related={'easix_day','siri_day','ldh_u_l','creatinine_mg_dl','platelets_109_l','wbc_109_l','neutrophils_harmonised','monocytes_harmonised','lymphocytes_harmonised',*meta.loc[meta.index_component,'variable'].tolist()}
features=base+meta.variable.tolist()
cohorts={}
for db in ['mimiciv','mimiciii']:
    d=pd.read_csv(W/'endpoints/data_private'/f'all_icu_d1d7_daily_{db}.csv',low_memory=False).query('day==1').copy()
    absolute=d[['neutrophils_abs','monocytes_abs','lymphocytes_abs']].notna().all(axis=1)&(db=='mimiciv')
    for c in ['neutrophils','monocytes','lymphocytes']:
        d[c+'_harmonised']=d[c+'_abs'].where(absolute,d.wbc_109_l*d[c+'_pct']/100)
    d['male']=d.sex.eq('M').astype(int)
    d=d.merge(pd.read_csv(W/'endpoints/data_private'/f'{db}_day1_system_means.csv'),on='icu_stay_id',validate='one_to_one')
    d=d.merge(labels.loc[labels.dataset.eq(db),['icu_stay_id','gbmtm_class']],on='icu_stay_id',validate='one_to_one')
    assert d.subject_id.is_unique
    cohorts[db]=d
iv=cohorts['mimiciv'];iii=cohorts['mimiciii']
train_idx,test_idx=train_test_split(np.arange(len(iv)),test_size=.30,random_state=SEED,stratify=iv.gbmtm_class)
train=iv.iloc[train_idx];internal=iv.iloc[test_idx]
assert set(train.subject_id).isdisjoint(internal.subject_id)
missing=train[features].isna().mean()
keep=[v for v in features if (v in related or missing[v]<.70) and train[v].nunique()>1]
x=train[keep];y=train.gbmtm_class.to_numpy()-1
availability=pd.DataFrame({'variable':features,'training_missing_pct':100*missing,'index_related':[v in related for v in features],'eligible':[v in keep for v in features]})
availability.to_csv(A/'training_feature_eligibility.csv',index=False)
parts=[]
for name,d in [('MIMIC-IV train',train),('MIMIC-IV internal validation',internal),('MIMIC-III external validation',iii)]:
    for c in [1,2,3]:parts.append(dict(cohort=name,trajectory_class=c,n=int(d.gbmtm_class.eq(c).sum())))
pd.DataFrame(parts).to_csv(A/'split_class_counts.csv',index=False)
pd.DataFrame({'icu_stay_id':iv.icu_stay_id,'subject_id':iv.subject_id,'gbmtm_class':iv.gbmtm_class,'partition':np.where(np.isin(np.arange(len(iv)),train_idx),'train','internal validation')}).to_csv(PRIVATE/'iv_patient_split.csv',index=False)

cv=StratifiedKFold(5,shuffle=True,random_state=SEED)
def prep(model):
    return Pipeline([('imputer',SimpleImputer(strategy='median')),
        ('quantile',QuantileTransformer(n_quantiles=100,output_distribution='normal',random_state=SEED)),('model',model)])

# 1. Sparse multinomial logistic model; refit each preprocessing step within CV.
l1=GridSearchCV(prep(LogisticRegression(solver='saga',l1_ratio=1.0,max_iter=10000,random_state=SEED)),
    {'model__C':[.01,.03,.1,.3,1,3]},scoring='roc_auc_ovr',cv=cv,n_jobs=3,refit=True)
l1.fit(x,y)
l1_coef=np.max(np.abs(l1.best_estimator_.named_steps['model'].coef_),axis=0)
s1=l1_coef>1e-6
print('L1 selected',int(s1.sum()),flush=True)

# 2. Shadow-feature RF: beat the strongest permuted feature in 40 independent fits.
z=prep(RandomForestClassifier(n_estimators=1)).steps[0][1].fit_transform(x)
shadow_hits=np.zeros(len(keep),int)
for iteration in range(40):
    shadows=np.column_stack([rng.permutation(z[:,j]) for j in range(z.shape[1])])
    rf=RandomForestClassifier(n_estimators=250,max_depth=6,class_weight='balanced',max_features='sqrt',random_state=SEED+iteration,n_jobs=3)
    rf.fit(np.column_stack([z,shadows]),y)
    importance=rf.feature_importances_
    shadow_hits+=importance[:len(keep)]>importance[len(keep):].max()
    if (iteration+1)%10==0:print('Shadow iterations',iteration+1,flush=True)
shadow_p=np.array([binomtest(int(h),40,.5,alternative='greater').pvalue for h in shadow_hits])
s2=shadow_p<.05/len(keep)
print('Shadow RF selected',int(s2.sum()),flush=True)

# 3. Recursive RF elimination, selected by five-fold train-only macro-AUROC.
rfe=RFECV(prep(RandomForestClassifier(n_estimators=200,max_depth=8,class_weight='balanced',random_state=SEED,n_jobs=1)),
    step=5,min_features_to_select=5,cv=cv,scoring='roc_auc_ovr',n_jobs=3,importance_getter='named_steps.model.feature_importances_')
rfe.fit(x,y);s3=rfe.support_
selected=[v for v,a,b,c in zip(keep,s1,s2,s3) if a and b and c]
selection=pd.DataFrame(dict(variable=keep,l1_coefficient=l1_coef,l1_selected=s1,shadow_hits=shadow_hits,shadow_p=shadow_p,shadow_selected=s2,rfe_rank=rfe.ranking_,rfe_selected=s3,consensus=[v in selected for v in keep]))
selection.to_csv(A/'ensemble_feature_selection.csv',index=False)
assert selected, 'Three-method intersection is empty; stop without substituting a union'
print('Consensus',len(selected),selected,flush=True)
pd.DataFrame(l1.cv_results_).drop(columns=['params']).to_csv(A/'l1_training_cv.csv',index=False)
pd.DataFrame({'n_features':rfe.cv_results_['n_features'],'mean_test_score':rfe.cv_results_['mean_test_score'],'std_test_score':rfe.cv_results_['std_test_score']}).to_csv(A/'rfe_training_cv.csv',index=False)

search=GridSearchCV(prep(XGBClassifier(objective='multi:softprob',num_class=3,tree_method='hist',eval_metric='mlogloss',learning_rate=.05,subsample=.85,colsample_bytree=.85,random_state=SEED,n_jobs=1)),
    {'model__max_depth':[2,3,4],'model__n_estimators':[200,400],'model__min_child_weight':[3,10]},cv=cv,scoring='roc_auc_ovr',n_jobs=3,refit=True)
search.fit(train[selected],y)
model=search.best_estimator_
pd.DataFrame(search.cv_results_).drop(columns=['params']).to_csv(A/'xgb_training_cv.csv',index=False)
joblib.dump(dict(features=selected,pipeline=model,seed=SEED),PRIVATE/'day1_classifier.joblib')
print('Fixed model fitted',search.best_params_,flush=True)

def metrics(yy,pp):
    onehot=np.eye(3)[yy]
    return dict(macro_auroc=roc_auc_score(yy,pp,multi_class='ovr'),macro_auprc=np.mean([average_precision_score(onehot[:,c],pp[:,c]) for c in range(3)]),
        accuracy=accuracy_score(yy,pp.argmax(axis=1)),balanced_accuracy=balanced_accuracy_score(yy,pp.argmax(axis=1)),macro_f1=f1_score(yy,pp.argmax(axis=1),average='macro'),
        log_loss=log_loss(yy,pp,labels=[0,1,2]),multiclass_brier=np.mean(np.sum((onehot-pp)**2,axis=1)))
perf=[];classperf=[];roc=[];cal=[];conf=[];shap_rows=[];predictions=[]
for name,d in [('MIMIC-IV train',train),('MIMIC-IV internal validation',internal),('MIMIC-III external validation',iii)]:
    yy=d.gbmtm_class.to_numpy()-1;pp=model.predict_proba(d[selected]);point=metrics(yy,pp)
    draws=[];perclass=[]
    # Stratified patient bootstrap retains all three classes in each replicate.
    for b in range(1000):
        ii=np.concatenate([rng.choice(np.flatnonzero(yy==c),sum(yy==c),replace=True) for c in range(3)])
        draws.append(metrics(yy[ii],pp[ii]))
        perclass.append([roc_auc_score((yy[ii]==c).astype(int),pp[ii,c]) for c in range(3)])
    samples=pd.DataFrame(draws)
    for m,value in point.items():
        lo,hi=samples[m].quantile([.025,.975]);perf.append(dict(cohort=name,n=len(d),metric=m,value=value,ci_low=lo,ci_high=hi,evaluation='apparent training fit' if name.endswith('train') else 'fixed-model held-out validation'))
    for c in range(3):
        lo,hi=np.quantile(np.array(perclass)[:,c],[.025,.975]);classperf.append(dict(cohort=name,trajectory_class=c+1,n=len(d),class_n=int(sum(yy==c)),auroc=roc_auc_score(yy==c,pp[:,c]),auroc_low=lo,auroc_high=hi,auprc=average_precision_score(yy==c,pp[:,c]),recall=float(np.mean(pp.argmax(axis=1)[yy==c]==c))))
        fpr,tpr,_=roc_curve(yy==c,pp[:,c])
        for f,t in zip(fpr,tpr):roc.append(dict(cohort=name,trajectory_class=c+1,fpr=f,tpr=t))
        bins=pd.qcut(pp[:,c],10,labels=False,duplicates='drop')
        for bn in np.unique(bins):
            mask=bins==bn;cal.append(dict(cohort=name,trajectory_class=c+1,bin=int(bn),n=int(mask.sum()),predicted=pp[mask,c].mean(),observed=np.mean(yy[mask]==c)))
    cc=confusion_matrix(yy,pp.argmax(axis=1),labels=[0,1,2])
    for i,j in np.ndindex(3,3):conf.append(dict(cohort=name,observed_class=i+1,predicted_class=j+1,n=int(cc[i,j]),proportion=cc[i,j]/cc[i].sum()))
    out=d[['icu_stay_id','subject_id','gbmtm_class']].copy();out['cohort']=name
    for c in range(3):out[f'prob_class_{c+1}']=pp[:,c]
    predictions.append(out)
    if not name.endswith('train'):
        tx=model[:-1].transform(d[selected])
        explainer=shap.TreeExplainer(model.named_steps['model'],feature_perturbation='tree_path_dependent',model_output='raw')
        explanation=explainer(tx)
        margin=model.named_steps['model'].predict(tx,output_margin=True)
        error=float(np.max(np.abs(explanation.base_values+explanation.values.sum(axis=1)-margin)))
        assert error<1e-4, error
        np.savez_compressed(PRIVATE/('shap_internal.npz' if 'internal' in name else 'shap_external.npz'),values=explanation.values,base_values=explanation.base_values,data=tx,features=selected)
        for j,v in enumerate(selected):
            for c in range(3):shap_rows.append(dict(cohort=name,variable=v,trajectory_class=c+1,n=len(d),mean_abs_shap=np.abs(explanation.values[:,j,c]).mean(),mean_shap=explanation.values[:,j,c].mean(),output_space='class-specific raw margin',max_additivity_error=error,background='tree-path-dependent training leaf counts'))
    print(name,point,flush=True)
for name,rows in [('prediction_performance',perf),('class_specific_performance',classperf),('roc_curves',roc),('calibration',cal),('confusion_matrix',conf),('shap_importance',shap_rows)]:
    pd.DataFrame(rows).to_csv(A/(name+'.csv'),index=False)
pd.concat(predictions).to_csv(PRIVATE/'day1_fixed_predictions.csv',index=False)
manifest=dict(seed=SEED,train_n=len(train),internal_validation_n=len(internal),external_validation_n=len(iii),patient_overlap_train_internal=0,
    candidate_n=len(features),eligible_n=len(keep),selector_counts={'l1':int(s1.sum()),'shadow_random_forest':int(s2.sum()),'rfecv':int(s3.sum())},selected_features=selected,
    design='MIMIC-IV patient-level stratified 70/30; fixed MIMIC-III historical cross-version external validation',
    l1_best_C=l1.best_params_,shadow_iterations=40,shadow_threshold='one-sided binomial p<0.05 / eligible feature count against strongest shadow',
    rfe_scoring='five-fold macro one-vs-rest AUROC within IV training only',best_xgboost_params=search.best_params_,train_cv_auroc=search.best_score_,
    preprocessing='train-fold median imputation, train-fold empirical quantile normalisation; final refit on IV training only',bootstrap='1000 stratified patient resamples, fixed model, percentile 95% intervals',
    versions=dict(sklearn=sklearn.__version__,xgboost=xgboost.__version__,shap=shap.__version__))
(A/'classifier_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Classifier, validation, bootstrap and SHAP complete.',flush=True)
