#!/usr/bin/env python3
"""Current article figures from traceable aggregates and saved probabilities."""
from pathlib import Path
import argparse, importlib.util, json, shutil, sys, tempfile
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np
import pandas as pd
from PIL import Image
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent
P = HERE.parents[1]
W = P/'src'
A = P/'data/aggregate/endpoints'
S = P/'data/aggregate/display'
FIXED = P/'data/aggregate/trajectory'
from figure_support import apply_style, export_figure, validate_flow, ENDPOINT_STYLES
apply_style()
COLORS = {1:'#4C78A8',2:'#59A14F',3:'#E15759'}
DB = {'mimiciv':'MIMIC-IV','mimiciii':'MIMIC-III'}
SYSTEM = {'Demography':'#8C939B','Severity':'#B697BF','Endothelial stress':'#D8A458',
          'Renal':'#74A8BC','Haematology':'#C78D9A','Inflammation':'#69A89A'}
GRAY, LIGHT = '#68717A', '#E6E8EB'
records=[]
OUTPUT = HERE
TREATMENT_SOURCE = S

def read(name): return pd.read_csv(S/(name+'.csv'))
def write(frame,name):
    path=S/(name+'.csv');frame.to_csv(path.with_suffix('.csv.tmp'),index=False)
    path.with_suffix('.csv.tmp').replace(path)
    return read(name)
def panel(ax,label):
    ax.text(-.12,1.09,label,transform=ax.transAxes,fontweight='bold',fontsize=10.5,va='top')
def clean(ax,grid='x'):
    ax.grid(axis=grid,color=LIGHT,lw=.5,zorder=0);ax.set_axisbelow(True)
def logaxis(ax,ticks,limits):
    ax.set_xscale('log');ax.set_xlim(*limits)
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.axvline(1,color=GRAY,lw=.7,ls='--');clean(ax)
def class_legend(fig):
    fig.legend(handles=[Line2D([0],[0],color=COLORS[c],ls=['-','--',':'][c-1],lw=1,label=f'Class {c}') for c in COLORS],
               loc='lower center',ncol=3,frameon=False)
def text_geometry(fig):
    """Check text-object bounds without opening or raster-inspecting the output."""
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    texts=list({id(t):t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()}.values())
    # Locators retain out-of-view tick objects; Axis.draw uses only _update_ticks().
    hidden_ticks=set()
    for ax in fig.axes:
        for axis in [ax.xaxis,ax.yaxis]:
            drawn={id(label) for tick in axis._update_ticks() for label in [tick.label1,tick.label2]}
            hidden_ticks.update(id(label) for tick in [*axis.get_major_ticks(),*axis.get_minor_ticks()]
                                for label in [tick.label1,tick.label2] if id(label) not in drawn)
    texts=[t for t in texts if id(t) not in hidden_ticks]
    bounds={id(t):t.get_window_extent(renderer) for t in texts}
    outside=[t.get_text() for t in texts if not (fig.bbox.contains(*bounds[id(t)].p0) and fig.bbox.contains(*bounds[id(t)].p1))]
    overlaps=[]
    for i,a in enumerate(texts):
        aa=bounds[id(a)]
        for b in texts[i+1:]:
            bb=bounds[id(b)]
            if min(aa.x1,bb.x1)-max(aa.x0,bb.x0)>.5 and min(aa.y1,bb.y1)-max(aa.y0,bb.y0)>.5:
                overlaps.append([a.get_text(),b.get_text()])
    for legend in [*fig.legends,*[ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]]:
        lb=legend.get_window_extent(renderer);own={id(t) for t in legend.get_texts()}
        for text in texts:
            if id(text) in own:continue
            tb=bounds[id(text)]
            if min(lb.x1,tb.x1)-max(lb.x0,tb.x0)>.5 and min(lb.y1,tb.y1)-max(lb.y0,tb.y0)>.5:
                overlaps.append(['legend area',text.get_text()])
    return dict(method='Matplotlib text/legend object bounds; no image QA',
        text_objects=len(texts),outside_canvas=outside,text_overlaps=overlaps)
def save(fig,stem):
    geometry=None
    if stem in {'figure_1_inclusion_exclusion_flow','figure_2_joint_trajectories','figure_3_class_treatment_outcomes','figure_4_class_clinical_profiles','figure_5_landmark_mortality_forest','figure_6_day7_landmark_survival','figure_7_day1_class_prediction','supp_figure_s2_day7_subgroups','supp_figure_s3_system_associations','supp_figure_s4_six_hour_management_profiles','supp_figure_s5_treatment_associations','supp_figure_s6_day1_management_profiles','supp_figure_s7_treatment_probability'}:
        geometry=text_geometry(fig)
        assert not geometry['outside_canvas'] and not geometry['text_overlaps'],geometry
    if geometry is None:geometry=text_geometry(fig)
    record=export_figure(fig,OUTPUT/stem)
    fig.savefig(OUTPUT/(stem+'.png'),dpi=300,facecolor='white',bbox_inches=None)
    record.update(stem=stem,class_colors=COLORS,axes_geometry=[dict(bounds=list(ax.get_position().bounds),title=ax.get_title(),xlabel=ax.get_xlabel(),ylabel=ax.get_ylabel()) for ax in fig.axes])
    if geometry:record['text_geometry']=geometry
    bounds=[a['bounds'] for a in record['axes_geometry']]
    record['panel_alignment']=dict(method='Numeric axes bounds; equal-sized panel groups',groups=[dict(width=k[0],height=k[1],panels=int(sum(round(b[2],6)==k[0] and round(b[3],6)==k[1] for b in bounds))) for k in sorted(set((round(b[2],6),round(b[3],6)) for b in bounds))])
    records.append(record);plt.close(fig)
def alias(source,target):
    for ext in ['.pdf','.tif','.png']:shutil.copy2(OUTPUT/(source+ext),OUTPUT/(target+ext))

def flow():
    path=HERE/'cohort_support.py'
    spec=importlib.util.spec_from_file_location('flow_skill',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    stem='figure_1_inclusion_exclusion_flow'
    cfg,steps=module.validate_config(json.loads((S/'figure_1_cohort_flow.json').read_text()))
    apply_style()
    width,height=183,225
    cfg['style']['height_mm']=height
    fig=plt.figure(figsize=(width/25.4,height/25.4));ax=fig.add_axes([0,0,1,1])
    ax.set(xlim=(0,width),ylim=(0,height));ax.axis('off')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    boxes=[]
    def box(x,y,w,h,text,fill='white',bold=False):
        wrapped=module.wrap_text(text,w-3,renderer,'Arial',8,bold)
        patch=module.FancyBboxPatch((x,y-h/2),w,h,boxstyle='round,pad=0,rounding_size=1.2',
            facecolor=fill,edgecolor='#202020',linewidth=.8)
        ax.add_patch(patch)
        ax.text(x+w/2,y,wrapped,ha='center',va='center',fontsize=8,
            fontweight='bold' if bold else 'normal',linespacing=1.3,parse_math=False)
        boxes.append(text)
    def arrow(a,b):
        ax.add_patch(module.FancyArrowPatch(a,b,arrowstyle='->',mutation_scale=8,
            color='#202020',linewidth=.8,shrinkA=0,shrinkB=0))
    for i,(co,group) in enumerate(zip(cfg['cohorts'],cfg['groups'])):
        left=2+i*91.5; mainw=42.5; sidew=41;side=left+mainw+4
        palette=module.PALETTES[group['palette']]
        box(left,202,87.5,10,group['label'],palette['header'],True)
        ax.text(left+43.75,190,co['label'],ha='center',va='center',fontweight='bold',fontsize=9)
        for j,node in enumerate(co['nodes']):
            y=173-29*j
            box(left,y,mainw,15,module.node_text(node),palette['fill'] if j==4 else 'white',j==4)
            if j<4:
                arrow((left+mainw/2,y-7.5),(left+mainw/2,y-29+7.5))
                by=y-14.5
                arrow((left+mainw/2,by),(side,by))
                box(side,by,sidew,12,module.node_text(co['exclusions'][j]),palette['fill'])
        # Analysis branches share the parent; availability is distinct from model N.
        analysis=co['analysis'];available=analysis['static_available'];models=analysis['static_models']
        sequence=lambda counts:' / '.join(f'{counts[key]:,}' for key in ['easix','siri','joint'])
        static=('Static · 24 h landmark\nEASIX / SIRI / Joint\nAvailable\n'+sequence(available)+
                '\n28-day mortality (Cox)\n'+sequence(models)+
                '\nICU mortality (FG)\n'+sequence(analysis['static_icu_models']))
        trajectory=(f'Trajectory analysis\nEligible / assigned: {analysis["trajectory_assigned"]:,}\n'
                    f'Day 7 landmark\n28-day mortality (Cox)\nN={analysis["trajectory_day7"]:,}\n'
                    f'ICU mortality (FG)\nN={analysis["trajectory_icu_day7"]:,}')
        box(left,22.5,mainw,40,static,palette['fill'])
        box(side,22.5,sidew,40,trajectory,palette['fill'])
        ax.plot([left+mainw/2,left+mainw/2,side+sidew/2],[49.5,44,44],color='#202020',lw=.8)
        arrow((left+mainw/2,44),(left+mainw/2,42.5))
        arrow((side+sidew/2,44),(side+sidew/2,42.5))
    save(fig,stem)
    audit=dict(cohorts=2,sequential_exclusions=len(steps),count_equations_passed=len(steps),
        zero_exclusions_shown=sum(s['excluded_n']==0 for s in steps),width_mm=width,height_mm=height,
        box_count=len(boxes),analysis_branches=4,bottom_explanatory_text=False,
        source='figure_source_data/figure_1_analysis_populations.csv',
        pdf_text='TrueType embedded',tif_dpi=300,tif_compression='tiff_lzw')
    records[-1].update(audit)
    (OUTPUT/(stem+'_source.json')).write_text(json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
    (OUTPUT/(stem+'_audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
    pd.DataFrame(steps).to_csv(OUTPUT/(stem+'_steps.csv'),index=False)

METRICS=[('nee_any','Any NEE','Patients (%)',True),('nee_high','NEE >0.3 µg/kg/min','Patients (%)',True),('fluid_kg_high','Net balance >30 mL/kg','Patients (%)',True),('mv','Mechanical ventilation','Patients (%)',True),('crrt','CRRT','Patients (%)',True),('map','Mean arterial pressure','mmHg',False),('urine','Urine output','mL/kg/h',False),('lactate','Blood lactate','mmol/L',False)]
def management():
    d=read('management_day1_summary')
    fig,axes=plt.subplots(4,2,figsize=(180/25.4,220/25.4))
    for i,(ax,(name,label,unit,binary)) in enumerate(zip(axes.flat,METRICS)):
        z=d[d.metric.eq(name)];mult=100 if binary else 1
        for db in DB:
            for r in z[z.dataset.eq(db)].itertuples():
                c=int(r.trajectory_class);y=3-c+(.12 if db=='mimiciv' else -.12)
                ax.errorbar(r.center*mult,y,xerr=[[(r.center-r.low)*mult],[(r.high-r.center)*mult]],fmt='o' if db=='mimiciv' else 'D',mfc=COLORS[c] if db=='mimiciv' else 'white',color=COLORS[c],capsize=2,ms=4,lw=1)
        ax.set(yticks=[2,1,0],yticklabels=['Class 1','Class 2','Class 3'],ylim=(-.6,2.6),xlabel=unit,xlim=(0,100) if binary else (0,None))
        ax.set_title(label,loc='left',fontweight='bold',fontsize=9)
        ps=[]
        for db in DB:
            pv=z[z.dataset.eq(db)].p_value.iloc[0];ps.append(DB[db]+(' P<0.001' if pv<.001 else f' P={pv:.3f}'))
        ax.text(0,.99,' · '.join(ps),transform=ax.transAxes,fontsize=7.5,va='top')
        clean(ax);panel(ax,chr(65+i))
    fig.legend(handles=[Line2D([0],[0],marker='o',color=GRAY,lw=0,label='MIMIC-IV'),Line2D([0],[0],marker='D',mfc='white',color=GRAY,lw=0,label='MIMIC-III')],loc='lower center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.13,right=.97,top=.94,bottom=.09,wspace=.48,hspace=.90)
    save(fig,'supp_figure_s6_day1_management_profiles')
    records[-1]['management']=dict(rows=len(d),metrics=8,source='figure_source_data/management_day1_summary.csv',binary_interval='Wilson 95% CI',continuous_interval='IQR')

def trajectories():
    tests=read('trajectory_daily_comparisons')
    fit=pd.read_csv(FIXED/'all_icu_gbmtm_selected_fitted_profiles.csv')
    fit['dataset']='mimiciv';fit['value']=fit.fitted_z;fit['statistic']='fitted model mean'
    external=pd.read_csv(A/'all_icu_gbmtm_trajectory_profiles.csv')
    external=external[external.dataset.eq('mimiciii')].copy()
    external['value']=external.mean_z;external['statistic']='mean preprocessed fixed-assignment profile'
    d=write(pd.concat([fit,external],ignore_index=True),'figure_2_selected_trajectories')
    counts=write(pd.read_csv(A/'all_icu_gbmtm_class_distribution.csv'),'trajectory_class_distribution')
    fig,axes=plt.subplots(2,2,figsize=(180/25.4,155/25.4))
    for i,(ax,(db,index)) in enumerate(zip(axes.flat,[(db,index) for db in DB for index in ['log2_easix_interp_z','log2_siri_interp_z']])):
        for c in COLORS:
            t=d[d.dataset.eq(db)&d['index'].eq(index)&d.gbmtm_class.eq(c)].sort_values('day')
            n=counts[counts.dataset.eq(db)&counts.gbmtm_class.eq(c)].iloc[0]
            ax.plot(t.day,t.value,color=COLORS[c],marker=['o','s','^'][c-1],ms=3.5,label=f'Class {c}: {int(n.n):,} ({n.pct:.1f}%)')
        ax.axhline(0,color=GRAY,lw=.6,ls=':')
        ax.set(xticks=range(1,8),xlim=(.65,7.35),ylim=(-1.6,1.95),xlabel='ICU day',ylabel='Standardised log2 index')
        name='EASIX' if 'easix' in index else 'SIRI'
        ax.set_title(f'{DB[db]} · {name}',loc='left',fontweight='bold')
        for r in tests[tests.dataset.eq(db)&tests['index'].eq(index)].itertuples():
            ax.text(r.day,1.70,r.stars,ha='center',va='center',fontsize=9.5,color='#202020')
        if name=='EASIX':
            handles,labels=ax.get_legend_handles_labels()
            fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.49 if db=='mimiciv' else .005),ncol=3,fontsize=8.5,handlelength=1,frameon=False,columnspacing=1.1)
        clean(ax,'y');panel(ax,chr(65+i))
    fig.subplots_adjust(left=.12,right=.98,top=.94,bottom=.15,hspace=.75,wspace=.40)
    save(fig,'figure_2_joint_trajectories')
    records[-1]['daily_comparisons']=dict(source='figure_source_data/trajectory_daily_comparisons.csv',
        days=list(range(1,8)),test='Kruskal-Wallis across the three fixed classes; observed values only',
        adjustment='BH over seven days within each database/index',
        stars={'*':'q<0.05','**':'q<0.01','***':'q<0.001','ns':'q>=0.05'})

SYSTEM_COLORS={'Cardiovascular':'#B98575','Respiratory':'#7DA7C4','Neurologic':'#A08CB6','Renal':'#74A8BC',
    'Hepatic':'#D1B071','Electrolytes':'#89AA84','Metabolic':'#BAA0A1','Haematology':'#C78D9A',
    'Coagulation':'#AC91A3','Index components':'#989EA7'}
SYSTEM_ORDER=list(SYSTEM_COLORS)

def phenotypes():
    d=read('baseline_comparison')
    metadata=d.drop_duplicates('variable').copy()
    metadata['order']=metadata.system.map({v:i for i,v in enumerate(SYSTEM_ORDER)})
    metadata=metadata.sort_values('order',kind='stable')
    variables=metadata.variable.tolist();n=len(variables)
    # Each class ranks its own IV means; III is paired to the same variable row.
    # Keep the system-ordered metadata above unchanged for supplementary Figure S3.
    fig=plt.figure(figsize=(180/25.4,250/25.4))
    axes=[fig.add_axes([.17+j/3,.38,.145,.52]) for j in range(3)]
    y=np.arange(n)[::-1];rank_rows=[]
    for ax,c in zip(axes,[1,2,3]):
        ranked=d[d.scope.eq('mimiciv')&d.state.eq(c)].sort_values('mean_z',ascending=False,kind='stable')
        assert len(ranked)==n and ranked.mean_z.is_monotonic_decreasing
        for j,r in enumerate(ranked.itertuples()):
            var=r.variable
            ext=d[d.scope.eq('mimiciii')&d.state.eq(c)&d.variable.eq(var)].iloc[0]
            ax.barh(y[j]+.10,r.mean_z,height=.46,color=COLORS[c],alpha=.82)
            ax.errorbar(r.mean_z,y[j]+.10,xerr=[[r.mean_z-r.ci_low],[r.ci_high-r.mean_z]],fmt='none',color='#29313B',lw=.5)
            ax.plot(ext.mean_z,y[j]-.23,'D',color=COLORS[c],mfc='white',ms=3,mew=.75)
            ax.axhspan(y[j]-.5,y[j]+.5,color=SYSTEM_COLORS[r.system],alpha=.07,zorder=-2)
            rank_rows.append(dict(trajectory_class=c,rank=j+1,variable=var,label=r.label,
                system=r.system,mimiciv_mean_z=r.mean_z,mimiciii_mean_z=ext.mean_z))
        ax.axvline(0,color=GRAY,lw=.6)
        ax.set(xlim=(-1.5,1.7),ylim=(-.7,n-.3),yticks=y,yticklabels=ranked.label,xticks=[-1,0,1],xlabel='Reference z score')
        ax.tick_params(axis='y',length=0,labelsize=8)
        ax.set_title(f'Class {c}',color=COLORS[c],fontweight='bold',pad=10);clean(ax)
        for tick,system in zip(ax.get_yticklabels(),ranked.system):tick.set_color(SYSTEM_COLORS[system])
    write(pd.DataFrame(rank_rows),'figure_4a_variable_order')
    fig.text(.02,.942,'A',fontweight='bold',fontsize=10.5)
    legend=[Patch(facecolor=v,label=k) for k,v in SYSTEM_COLORS.items() if k in metadata.system.tolist()]
    fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.999),ncol=5,fontsize=8.2,handlelength=.8,columnspacing=.8,frameon=False)
    fig.legend(handles=[Patch(facecolor=GRAY,label='MIMIC-IV'),
        Line2D([0],[0],marker='D',markerfacecolor='white',markeredgecolor=GRAY,
               color='none',markersize=4,label='MIMIC-III')],
        loc='upper center',bbox_to_anchor=(.63,.953),ncol=2,fontsize=8.5,frameon=False)
    radar=read('subtype_system_radar')
    systems=[s for s in SYSTEM_ORDER if s!='Index components']
    angles=np.linspace(0,2*np.pi,len(systems),endpoint=False)
    bottom=fig.add_gridspec(1,2,left=.12,right=.90,top=.275,bottom=.07,wspace=.65)
    for i,db in enumerate(DB):
        ax=fig.add_subplot(bottom[0,i],projection='polar');ax.set_theta_offset(np.pi/2);ax.set_theta_direction(-1)
        for c in COLORS:
            values=radar[radar.scope.eq(db)&radar.state.eq(c)].set_index('system').loc[systems].mean_percentile.to_numpy()
            assert np.all((values>=30)&(values<=70))
            ax.plot(np.r_[angles,angles[0]],np.r_[values,values[0]],color=COLORS[c],
                marker=['o','s','^'][c-1],ms=3,lw=1.4)
        ax.set_xticks(angles,systems,fontsize=8.5)
        ax.set(ylim=(30,70),yticks=[40,50,60,70],yticklabels=['40','50','60','70'])
        ax.tick_params(axis='y',labelsize=7.5);ax.set_rlabel_position(10)
        ax.plot(np.linspace(0,2*np.pi,181),np.full(181,50),ls=':',color=GRAY,lw=.75,zorder=0)
        ax.grid(color=LIGHT,lw=.6);ax.spines['polar'].set_color(LIGHT)
        ax.set_title(DB[db],fontsize=9,fontweight='bold',pad=16)
    fig.text(.02,.30,'B',fontweight='bold',fontsize=10.5)
    fig.legend(handles=[Line2D([0],[0],color=COLORS[c],marker=['o','s','^'][c-1],lw=1.4,label=f'Class {c}') for c in COLORS],loc='lower center',bbox_to_anchor=(.5,.001),ncol=3,fontsize=9,frameon=False)
    save(fig,'figure_4_class_clinical_profiles')
    records[-1]['cohort_legend']={'MIMIC-IV':'solid bar','MIMIC-III':'open diamond',
        'placement':'Dedicated panel A header band above the class panels'}
    records[-1]['variable_order']=dict(rule='Descending MIMIC-IV observed mean z score within each class',
        external_alignment='MIMIC-III follows each class-specific MIMIC-IV variable order',
        source='figure_source_data/figure_4a_variable_order.csv',rows=len(rank_rows))
    records[-1]['radar']=dict(axes=systems,datasets=list(DB),classes=[1,2,3],scale=[30,70],
        reference_midpoint=50,display='Two enlarged database panels overlaying all three classes',
        source='figure_source_data/subtype_system_radar.csv',
        mapping='figure_source_data/system_radar_variables.csv',
        aggregation='Equal mean of observed variable mean MIMIC-IV reference percentiles',
        direction='Measured-value direction retained; no severity reorientation')
    # Correlation display retains all selected systems, including external overlap.
    corr=read('system_associations');fig,axes=plt.subplots(1,2,figsize=(180/25.4,205/25.4))
    for i,(ax,db) in enumerate(zip(axes,DB)):
        matrix=corr[corr.scope.eq(db)].pivot(index='variable',columns='state',values='r').loc[variables,[1,2,3]].to_numpy()
        im=ax.imshow(matrix,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
        for yy,xx in np.ndindex(matrix.shape):
            val=matrix[yy,xx];ax.text(xx,yy,f'{val:+.2f}',ha='center',va='center',fontsize=7.5,color='white' if abs(val)>.65 else '#202020')
        ax.set_xticks([0,1,2],['Class 1','Class 2','Class 3'])
        for tick,c in zip(ax.get_xticklabels(),COLORS):tick.set_color(COLORS[c])
        ax.set_yticks(range(n),metadata.label if i==0 else ['']*n);ax.tick_params(length=0,labelsize=8)
        ax.set_title(DB[db],loc='left',fontweight='bold')
        ax.text(-.12,1.04,chr(65+i),transform=ax.transAxes,fontweight='bold',fontsize=10.5,va='top')
        for spine in ax.spines.values():spine.set_visible(False)
    fig.subplots_adjust(left=.27,right=.87,top=.94,bottom=.08,wspace=.16)
    colorbar=fig.colorbar(im,cax=fig.add_axes([.89,.18,.015,.65]),label='Point-biserial r')
    colorbar.ax.yaxis.label.set_fontsize(8.5);colorbar.ax.tick_params(labelsize=8)
    save(fig,'supp_figure_s3_system_associations')

def mortality_forest():
    d=read('figure_5_mortality_forest').copy()
    selections=[('static_joint','log2_easix_z','EASIX'),('static_joint','log2_siri_z','SIRI'),
        ('gbmtm_k3_d2_pruning','trajectory_class2','Class 2 vs 1'),('gbmtm_k3_d2_pruning','trajectory_class3','Class 3 vs 1')]
    d['row']=d.term.map({v[1]:i for i,v in enumerate(selections)})+d.dataset.map({'mimiciv':0,'mimiciii':6})
    assert len(d)==16 and set(d.horizon)=={'28-day','ICU'}
    fig=plt.figure(figsize=(180/25.4,135/25.4))
    axes=[fig.add_axes([.22,.13,.14,.73]),fig.add_axes([.63,.13,.14,.73])]
    for i,(ax,h) in enumerate(zip(axes,['28-day','ICU'])):
        for r in d[d.horizon.eq(h)].itertuples():
            color=COLORS[int(r.term[-1])] if r.term.startswith('trajectory') else '#526678';yy=10-r.row
            ax.errorbar(r.hazard_ratio,yy,xerr=[[r.hazard_ratio-r.ci_low],[r.ci_high-r.hazard_ratio]],fmt='s' if r.term.startswith('trajectory') else 'o',color=color,mfc=color if r.dataset=='mimiciv' else 'white',ms=4,capsize=2,lw=1)
            pv='<0.001' if r.p_value<.001 else f'={r.p_value:.3f}'
            ax.text(1.08,yy,f'{r.hazard_ratio:.2f} ({r.ci_low:.2f}–{r.ci_high:.2f})\nP{pv}',transform=ax.get_yaxis_transform(),fontsize=8,va='center',linespacing=1.15)
        logaxis(ax,[1,2,3],(.95,3.1));ax.set(ylim=(-.3,11.8),yticks=[],xlabel='Adjusted sHR' if h=='ICU' else 'Adjusted HR')
        ax.set_title('ICU mortality' if h=='ICU' else '28-day mortality',fontweight='bold',pad=15);ax.text(-.18,1.14,chr(65+i),transform=ax.transAxes,fontweight='bold',fontsize=10.5)
    for db,offset in [('mimiciv',0),('mimiciii',6)]:
        axes[0].text(-.12,11.1-offset,DB[db],ha='right',transform=axes[0].get_yaxis_transform(),fontweight='bold',fontsize=8.5)
        for j,(_,term,label) in enumerate(selections):
            axes[0].text(-.12,10-j-offset,label,ha='right',va='center',transform=axes[0].get_yaxis_transform(),fontsize=8.5,color=COLORS[int(term[-1])] if term.startswith('trajectory') else '#202020')
    save(fig,'figure_5_landmark_mortality_forest')

def km():
    curves,risk,tests=[pd.read_csv(A/n) for n in ['all_icu_km_gbmtm_curve_data.csv','all_icu_km_gbmtm_numbers_at_risk.csv','all_icu_km_gbmtm_logrank_tests.csv']]
    curves,risk,tests=[x[x.analysis.eq('day7_landmark')&x.horizon.eq('28-day')].copy() for x in [curves,risk,tests]]
    base=P/'data/aggregate/icu'
    cif,cr,ct=[pd.read_csv(base/n) for n in ['cif_display_curve_data.csv','cif_display_numbers_at_risk.csv','cif_display_gray_tests.csv']]
    assert set(cif.event)=={1}
    fig=plt.figure(figsize=(180/25.4,190/25.4));outer=fig.add_gridspec(2,2,hspace=.72,wspace=.42)
    alignment=[]
    for i,(db,h) in enumerate([(db,h) for db in DB for h in ['28-day','ICU']]):
        grid=outer[i//2,i%2].subgridspec(2,1,height_ratios=[3,1.1],hspace=.40)
        ax,ra=fig.add_subplot(grid[0]),fig.add_subplot(grid[1]);ra.axis('off')
        ra.text(0,1.12,'Number at risk',fontsize=8.5,transform=ra.transAxes)
        ticks=[7,14,21,28] if h=='28-day' else sorted(cr.admission_day.unique())
        for c in COLORS:
            t=(curves if h=='28-day' else cif);t=t[t.dataset.eq(db)&t.trajectory_class.eq(c)]
            x=t.time+7 if h=='28-day' else t.admission_day
            y=t.survival if h=='28-day' else t.cif
            ax.step(x,y,where='post',color=COLORS[c],lw=1.2,ls=['-','--',':'][c-1])
            ax.fill_between(x,t.ci_low,t.ci_high,step='post',color=COLORS[c],alpha=.1,lw=0)
            rr=(risk if h=='28-day' else cr);rr=rr[rr.dataset.eq(db)&rr.trajectory_class.eq(c)].copy()
            if h=='28-day':rr['admission_day']=rr.time+7;rr['event_free_icu_risk_n']=rr.n_risk
            rr=rr.set_index('admission_day')
            ra.text(-.06,.82-(c-1)*.3,f'C{c}',color=COLORS[c],ha='right',fontsize=8.5,transform=ra.transAxes)
            for tt in ticks:ra.text((tt-7)/(ticks[-1]-7),.82-(c-1)*.3,f'{int(rr.loc[tt,"event_free_icu_risk_n"]):,}',ha='center',fontsize=8.5,transform=ra.transAxes)
        test=(tests if h=='28-day' else ct);test=test[test.dataset.eq(db)].iloc[0]
        pv=test.logrank_p if h=='28-day' else test.p_value
        ax.text(.98,.05,('Log-rank P' if h=='28-day' else 'Gray P')+('<0.001' if pv<.001 else f'={pv:.3f}'),ha='right',transform=ax.transAxes,fontsize=8.5)
        ax.set(xlim=(7,ticks[-1]),ylim=(.3,1.01) if h=='28-day' else (0,.65),xticks=ticks,ylabel='Survival probability' if h=='28-day' else 'Cumulative incidence')
        ra.text(.5,-.34,'Days since ICU admission',ha='center',fontsize=8.5,transform=ra.transAxes)
        ax.set_title(f'{DB[db]}\n'+('28-day mortality' if h=='28-day' else 'ICU mortality'),loc='left',fontweight='bold');clean(ax,'y');ax.text(-.17,1.17,chr(65+i),transform=ax.transAxes,fontweight='bold',fontsize=10.5)
        alignment.append(dict(panel=chr(65+i),curve_bounds=list(ax.get_position().bounds),risk_bounds=list(ra.get_position().bounds)))
    class_legend(fig);fig.subplots_adjust(left=.12,right=.96,top=.92,bottom=.14)
    save(fig,'figure_6_day7_landmark_survival')
    records[-1]['survival']=dict(km_horizon='28-day',cif_days=[7,60],cif_event=1,competing_event='Live ICU discharge',confidence_interval=.95,alignment=alignment)

def subgroup():
    d=read('figure_7_subgroup_effects');q=read('figure_7_subgroup_interactions')
    order=['sepsis','aki','age','sex']
    strata=[['No sepsis','Sepsis'],['No AKI','AKI'],['Age <65','Age >=65'],['Female','Male']]
    fig=plt.figure(figsize=(180/25.4,190/25.4))
    for i,(db,h) in enumerate([(db,h) for db in DB for h in ['28-day','ICU']]):
        ax=fig.add_axes([.17+.49*(i%2),.59 if i<2 else .12,.15,.31]);labels=[v.replace('>=','≥') for ss in strata for v in ss]
        for j,(group,ss) in enumerate(zip(order,strata)):
            for k,st in enumerate(ss):
                y=7-j*2-k;z=d[d.dataset.eq(db)&d.horizon.eq(h)&d.subgroup.eq(group)&d.stratum.eq(st)]
                for c in [2,3]:
                    t=z[z.term.eq(f'trajectory_class{c}')]
                    if t.empty or not np.isfinite(t.iloc[0].hazard_ratio):continue
                    r=t.iloc[0];ax.errorbar(r.hazard_ratio,y+(.13 if c==2 else -.13),xerr=[[r.hazard_ratio-r.ci_low],[r.ci_high-r.hazard_ratio]],color=COLORS[c],fmt='o' if c==2 else 's',ms=3.8,capsize=1.8,lw=.9)
                label='Not estimable' if z.empty else f'{int(z.iloc[0].n):,}/{int(z.iloc[0].events):,}'
                ax.text(1.03,y,label,transform=ax.get_yaxis_transform(),va='center',fontsize=8.5)
            qq=q[q.dataset.eq(db)&q.horizon.eq(h)&q.subgroup.eq(group)].iloc[0]
            ax.text(1.65,6.5-j*2,f'{qq.interaction_p_bh:.3f}',transform=ax.get_yaxis_transform(),va='center',fontsize=7.5)
        ax.set(yticks=range(8),yticklabels=labels[::-1],ylim=(-.6,7.6),xlabel=('Adjusted sHR' if h=='ICU' else 'Adjusted HR')+' (95% CI)')
        logaxis(ax,[.5,1,2,4,6],(.45,6.5))
        ax.text(1.03,1.04,'N/events',transform=ax.transAxes,fontsize=8.5,fontweight='bold')
        ax.text(1.65,1.04,'BH P',transform=ax.transAxes,fontsize=8.5,fontweight='bold')
        ax.tick_params(axis='y',labelsize=7.5)
        ax.set_title(DB[db]+'\n'+('ICU mortality' if h=='ICU' else '28-day mortality'),loc='left',fontweight='bold',pad=29,fontsize=9)
        ax.text(-.12,1.22,chr(65+i),transform=ax.transAxes,fontweight='bold',fontsize=10.5)
    fig.legend(handles=[Line2D([0],[0],marker='o',color=COLORS[2],lw=0,label='Class 2 vs 1'),Line2D([0],[0],marker='s',color=COLORS[3],lw=0,label='Class 3 vs 1')],loc='lower center',ncol=2,frameon=False)
    save(fig,'supp_figure_s2_day7_subgroups')
    records[-1]['subgroups']=dict(displayed=order,omitted_treatments=['ventilation','vasopressor','crrt'],
        bh_family_tests=7,correction='Unchanged seven-test family separately by database and outcome',models_refitted=False)

def classifier():
    import shap
    source=P/'data/aggregate/classifier'
    perf=pd.read_csv(source/'prediction_performance.csv');roc=pd.read_csv(source/'roc_curves.csv')
    cal=pd.read_csv(source/'calibration.csv');shapdata=pd.read_csv(source/'shap_importance.csv')
    for data,name in [(perf,'performance'),(roc,'roc'),(cal,'calibration'),(shapdata,'shap')]:write(data,'figure_7_prediction_'+name)
    cohorts=['MIMIC-IV internal validation','MIMIC-III external validation']
    fig=plt.figure(figsize=(180/25.4,235/25.4))
    # Validation reads down the left column; one full beeswarm per output on the right.
    validation_axes=[fig.add_axes([.12,y,.25,.145]) for y in [.795,.565,.335,.105]]
    for i,cohort in enumerate(cohorts):
        ax=validation_axes[i];ax.plot([0,1],[0,1],ls=':',color=GRAY,lw=.7)
        for c in COLORS:
            t=roc[roc.cohort.eq(cohort)&roc.trajectory_class.eq(c)]
            ax.plot(t.fpr,t.tpr,color=COLORS[c],lw=1,label=f'Class {c}')
        val=perf[perf.cohort.eq(cohort)&perf.metric.eq('macro_auroc')].iloc[0]
        ax.text(.96,.07,f'Macro-AUROC {val.value:.3f}\n95% CI {val.ci_low:.3f}–{val.ci_high:.3f}',ha='right',transform=ax.transAxes,fontsize=9.3)
        ax.set(xlim=(0,1),ylim=(0,1),xlabel='False-positive rate',ylabel='True-positive rate',xticks=[0,.5,1],yticks=[0,.5,1])
        ax.set_title('MIMIC-IV internal · ROC' if i==0 else 'MIMIC-III external · ROC',loc='left',fontweight='bold',fontsize=10);clean(ax,'both')
    ax=validation_axes[2];metrics=['macro_auroc','macro_auprc','balanced_accuracy']
    for i,(cohort,color) in enumerate(zip(cohorts,['#526678','#D8A458'])):
        t=perf[perf.cohort.eq(cohort)].set_index('metric').loc[metrics]
        ax.bar(np.arange(3)+(i-.5)*.34,t.value,width=.31,color=color,label='MIMIC-IV internal' if i==0 else 'MIMIC-III external')
        ax.errorbar(np.arange(3)+(i-.5)*.34,t.value,yerr=[t.value-t.ci_low,t.ci_high-t.value],color='#202020',fmt='none',capsize=2,lw=.7)
    ax.set(ylim=(0,1.45),xticks=np.arange(3),xticklabels=['Macro-\nAUROC','Macro-\nAUPRC','Balanced\naccuracy'],ylabel='Score',yticks=[0,.5,1])
    ax.set_title('Validation scores',loc='left',fontweight='bold',fontsize=10)
    score_legend=ax.legend(loc='upper right',fontsize=8.7,frameon=False)
    fig.canvas.draw()
    highest_ci=perf[perf.cohort.isin(cohorts)&perf.metric.isin(metrics)].ci_high.max()
    clearance=(score_legend.get_window_extent(fig.canvas.get_renderer()).y0-
               ax.transData.transform((0,highest_ci))[1])*72/fig.dpi
    assert clearance>3, f'Figure 7C legend requires more headroom: {clearance:.2f} pt'
    score_legend_info=dict(location='upper right',reserved_headroom=True,
        minimum_vertical_clearance_pt=float(clearance),data_score_ticks=[0,.5,1])
    clean(ax,'y')
    ax=validation_axes[3];ax.plot([0,1],[0,1],ls=':',color=GRAY,lw=.7)
    for c in COLORS:
        for cohort,ls,marker in zip(cohorts,['-','--'],['o','s']):
            t=cal[cal.cohort.eq(cohort)&cal.trajectory_class.eq(c)]
            ax.plot(t.predicted,t.observed,color=COLORS[c],ls=ls,marker=marker,ms=2.5,lw=.9)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='Predicted probability',ylabel='Observed proportion',xticks=[0,.5,1],yticks=[0,.5,1])
    ax.set_title('Calibration',loc='left',fontweight='bold',fontsize=10);clean(ax,'both')
    for i,ax in enumerate(validation_axes):
        ax.tick_params(axis='both',labelsize=9.5)
        ax.xaxis.label.set_fontsize(10);ax.yaxis.label.set_fontsize(10)
        fig.text(.025,ax.get_position().y1+.017,chr(65+i),fontsize=10.5,fontweight='bold')
    labels=read('display_variables').set_index('variable').label.to_dict()
    labels.update(easix_day='EASIX',siri_day='SIRI',platelets_109_l='Platelets · minimum',
                  monocytes_harmonised='Monocytes · index',platelet_mean='Platelets · mean',baseexcess_mean='Base excess')
    with np.load(W/'classifier/data_private/shap_internal.npz',allow_pickle=False) as saved:
        features=saved['features'].tolist()
        values=saved['values'].copy();base_values=saved['base_values'].copy();data=saved['data'].copy()
    manifest=json.loads((source/'classifier_manifest.json').read_text())
    assert features==manifest['selected_features']
    assert values.shape==(manifest['internal_validation_n'],len(features),3)
    explanation=shap.Explanation(values=values,base_values=base_values,data=data,
        feature_names=[labels[v] for v in features],output_names=['Class 1','Class 2','Class 3'])
    rows=[];panel_records=[]
    np.random.seed(20261111)
    for c,y in zip(COLORS,[.72,.415,.11]):
        ax=fig.add_axes([.68,y,.29,.22])
        class_exp=explanation[...,c-1]
        # Match the saved float32 column-wise aggregation, avoiding layout-dependent summation rounding.
        importance=np.array([np.abs(values[:,j,c-1]).mean() for j in range(len(features))])
        expected=shapdata[shapdata.cohort.eq(cohorts[0])&shapdata.trajectory_class.eq(c)].set_index('variable').loc[features]
        np.testing.assert_allclose(importance,expected.mean_abs_shap,rtol=1e-6,atol=1e-8)
        order=np.argsort(importance)[::-1][:8]
        shap.plots.beeswarm(class_exp,max_display=8,order=shap.Explanation.abs.mean(0),
            ax=ax,show=False,plot_size=None,group_remaining_features=False,
            s=9,alpha=.8,color_bar=True,color_bar_label='Feature value')
        ax.set_title(f'SHAP · Class {c}',loc='left',color=COLORS[c],fontweight='bold',fontsize=10,pad=10)
        ax.set_xlabel('SHAP value (raw margin)',fontsize=10)
        ax.tick_params(axis='both',labelsize=10.4)
        letter=chr(68+c)
        fig.text(.47,y+.25,letter,fontweight='bold',fontsize=10.5)
        for rank,j in enumerate(order,1):
            rows.append(dict(panel=letter,trajectory_class=c,rank=rank,variable=features[j],
                label=labels[features[j]],mean_abs_shap=float(importance[j]),n=len(values),
                cohort=cohorts[0],output_space='class-specific raw margin',
                selection='all internal-validation patients for each model output'))
        panel_records.append(dict(panel=letter,trajectory_class=c,output_index=c-1,n=len(values),
            displayed_features=[features[j] for j in order],class_color=COLORS[c]))
    # The library colorbars keep SHAP's standard Low/High encoding, independent of class color.
    for ax in fig.axes:
        if ax.get_ylabel()=='Feature value':
            ax.yaxis.label.set_fontsize(10);ax.tick_params(labelsize=9.5)
    write(pd.DataFrame(rows),'figure_7_shap_panels')
    fig.legend(handles=[Line2D([0],[0],color=COLORS[c],lw=2,label=f'Class {c}') for c in COLORS]+[Line2D([0],[0],color=GRAY,ls='-',marker='o',label='MIMIC-IV internal'),Line2D([0],[0],color=GRAY,ls='--',marker='s',label='MIMIC-III external')],loc='lower center',ncol=3,fontsize=9.8,frameon=False)
    save(fig,'figure_7_day1_class_prediction')
    records[-1]['score_legend']=score_legend_info
    records[-1]['shap']=dict(plot_api='shap.plots.beeswarm',shap_version=shap.__version__,
        seed=20261111,panels=panel_records,features_per_panel=8,group_remaining_features=False,
        ranking='Class-specific mean absolute SHAP value across all 810 MIMIC-IV internal-validation patients',
        feature_color='Standard SHAP red-blue: high-low transformed model feature value',
        background='tree-path-dependent training leaf counts',output_space='class-specific raw margin',
        aggregate_source='figure_source_data/figure_7_shap_panels.csv')

def supporting():
    d=pd.read_csv(FIXED/'all_icu_gbmtm_candidate_metrics_k1_k6.csv');write(d,'supp_figure_s1_k_selection')
    fig,axes=plt.subplots(2,2,figsize=(180/25.4,130/25.4))
    for i,(ax,(field,label,threshold)) in enumerate(zip(axes.flat,[('BIC','BIC',None),('entropy','Entropy',.8),('min_avepp','Minimum AvePP',.7),('min_class_pct','Smallest class (%)',10)])):
        ax.plot(d.k,d[field],color='#526678',marker='o',ms=4);ax.axvline(3,color='#D8A458',ls='--',lw=1)
        if threshold is not None:ax.axhline(threshold,color=GRAY,ls=':',lw=.8)
        ax.set(xticks=range(1,7),xlabel='Number of classes (K)',ylabel=label);clean(ax,'y');panel(ax,chr(65+i))
    fig.subplots_adjust(left=.12,right=.98,top=.92,bottom=.12,wspace=.40,hspace=.50);save(fig,'supp_figure_s1_gbmtm_k_selection')
def management_profiles():
    d=read('supp_management_profiles')
    metrics=[(n,l+'\n('+u+')',100 if b else 1) for n,l,u,b in METRICS]
    fig,axes=plt.subplots(8,3,figsize=(180/25.4,250/25.4),sharex=True)
    for row,(metric,label,mult) in enumerate(metrics):
        z=d[d.metric.eq(metric)]
        low=z.ci_lower if metric in ['nee_any','nee_high','fluid_kg_high','mv','crrt'] else z.p25;high=z.ci_upper if metric in ['nee_any','nee_high','fluid_kg_high','mv','crrt'] else z.p75
        ymin=low.min()*mult;ymax=high.max()*mult
        pad=(ymax-ymin)*.12
        for c in COLORS:
            ax=axes[row,c-1]
            for db in DB:
                t=z[z.trajectory_class.eq(c)&z.dataset.eq(db)].sort_values('time_start_h')
                lo=t.ci_lower if metric in ['nee_any','nee_high','fluid_kg_high','mv','crrt'] else t.p25;hi=t.ci_upper if metric in ['nee_any','nee_high','fluid_kg_high','mv','crrt'] else t.p75
                ax.errorbar(np.arange(4)+(-.035 if db=='mimiciv' else .035),t.center*mult,
                    yerr=[np.maximum(0,(t.center-lo)*mult),np.maximum(0,(hi-t.center)*mult)],color=COLORS[c],
                    marker='o' if db=='mimiciv' else 'D',linestyle='-' if db=='mimiciv' else '--',
                    markerfacecolor=COLORS[c] if db=='mimiciv' else 'white',ms=3.2,capsize=1.8,lw=.9)
            ax.set(ylim=(max(0,ymin-pad),ymax+pad),xlim=(-.25,3.25),xticks=range(4),xticklabels=['0–6','6–12','12–18','18–24'])
            if c==1:ax.set_ylabel(label.replace('NEE >0.3 µg/kg/min','NEE >0.3\nµg/kg/min').replace('Net balance >30 mL/kg','Net balance\n>30 mL/kg').replace('Mechanical ventilation','Mechanical\nventilation'),rotation=0,ha='right',va='center',fontsize=8.5)
            if row==0:ax.set_title(f'Class {c}',fontweight='bold',color=COLORS[c],pad=10)
            if row==7:ax.set_xlabel('ICU hours')
            clean(ax,'y')
    fig.legend(handles=[Line2D([0],[0],color=GRAY,marker='o',lw=.9,label='MIMIC-IV'),
        Line2D([0],[0],color=GRAY,marker='D',mfc='white',ls='--',lw=.9,label='MIMIC-III')],
        loc='lower center',bbox_to_anchor=(.57,.002),ncol=2,frameon=False,fontsize=8.5)
    fig.subplots_adjust(left=.27,right=.98,top=.94,bottom=.09,hspace=.36,wspace=.23)
    save(fig,'supp_figure_s4_six_hour_management_profiles')
    records[-1]['management_profiles']=dict(datasets=list(DB),metrics=8,rows=len(d),
        class_panels=3,database_markers={'mimiciv':'solid circle','mimiciii':'open diamond; dashed line'},
        source='figure_source_data/supp_management_profiles.csv')

def treatment_forest():
    from treatment_figures import fixed_class_forest
    save(fixed_class_forest(TREATMENT_SOURCE), 'supp_figure_s5_treatment_associations')
    source = TREATMENT_SOURCE / 'figure_s6_fixed.csv'
    estimates = pd.read_csv(source)
    estimated = int(estimates.status.eq('estimated').sum())
    records[-1]['treatment_associations'] = dict(
        estimates=len(estimates), estimated=estimated, not_estimable=len(estimates)-estimated,
        global_tests=len(estimates.groupby(['dataset','outcome','exposure'])),
        datasets=list(DB), horizons=['28-day','ICU'], classes=[1,2,3],
        x_scale='log', x_limits=[.2,6], effect_source=str(source.relative_to(P)),
        interaction_source=str(source.relative_to(P)), nominal_p_values=True, models_refitted=False,
        exposures=list(estimates.exposure.unique()), encoding='clinical cutoff binary',
        truncation='Arrowheads indicate intervals outside displayed limits; text reports untruncated estimates')


def treatment_probability():
    from treatment_figures import treatment_probability as build_figure
    save(build_figure(TREATMENT_SOURCE), 'supp_figure_s7_treatment_probability')
    records[-1]['treatment_probability'] = dict(
        sources=[str((TREATMENT_SOURCE / f'figure_s7_four_treatment_{part}.csv').relative_to(P))
                 for part in ('curves','tests','histograms')],
        nominal_p_values=True, models_refitted=False)


def exact_alluvial(ax,nodes,links,stages,stage_categories,styles,label_mode='all',bar_width=None,small_label_threshold=0):
    """Draw observed counts, preserving source-class subflows through every node."""
    categories=list(dict.fromkeys(category for group in stage_categories for category in group))
    check=validate_flow(nodes,links,stages,categories=categories)
    n=check['cohort_stays'];gap=.035*n
    height=n+(max(map(len,stage_categories))-1)*gap
    positions={};bar=bar_width if bar_width is not None else (.29 if label_mode=='counts' else .64)
    ax.set(xlim=(-.36,len(stages)-.64),ylim=(height,0));ax.set_axis_off()
    renderer=ax.figure.canvas.get_renderer();annotations=[];outside_labels=[]
    for x,(stage,group) in enumerate(zip(stages,stage_categories)):
        source_nodes=nodes[nodes.hour.eq(stage)].set_index('category')
        counts=source_nodes.stays
        active=[key for key in group if counts.get(key,0)>0]
        cursor=(height-n-max(0,len(active)-1)*gap)/2
        for category in active:
            count=int(counts[category]);positions[stage,category]=(cursor,count)
            style=styles[category]
            patch=Rectangle((x-bar/2,cursor),bar,count,facecolor=style['color'],
                edgecolor=style.get('edgecolor',style['color']),hatch=style.get('hatch',''),linewidth=.45,zorder=4)
            ax.add_patch(patch)
            row=source_nodes.loc[category]
            assert np.isclose(row.panel_percent,100*count/n) and row.percentage_label==f'{100*count/n:.1f}%'
            rgb=np.array(to_rgb(style['color']))
            linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
            color='white' if linear@np.array([.2126,.7152,.0722])<.179 else '#111111'
            label=row.percentage_label if label_mode=='counts' else f'{style.get("label",category)}\n{count:,}\n{row.percentage_label}'
            text=ax.text(x,cursor+count/2,label,ha='center',va='center',fontsize=8.5,
                linespacing=1.0,color=color,zorder=6)
            bounds=text.get_window_extent(renderer);node_bounds=patch.get_window_extent(renderer)
            inside=bounds.width+2<=node_bounds.width and bounds.height+2<=node_bounds.height
            if label_mode=='counts' and row.panel_percent<small_label_threshold:inside=False
            placement='inside'
            if label_mode=='counts' and not inside:
                # Preserve true tiny-node heights; place their percentage beside the node.
                left=category=='unknown' and small_label_threshold==0;offset=bar/2+.05
                text.set_position((x-offset if left else x+offset,cursor+count/2))
                text.set_ha('right' if left else 'left');text.set_color('#202020')
                text.set_bbox(dict(facecolor='white',edgecolor='none',pad=.2))
                ax.plot([x-bar/2 if left else x+bar/2,x-offset+.012 if left else x+offset-.012],
                    [cursor+count/2]*2,color='#434A50',lw=.5,zorder=5)
                placement='adjacent to narrow segment'
                outside_labels.append((x,text,cursor+count/2))
            else:
                assert inside,(category,label,bounds.width,bounds.height,node_bounds.width,node_bounds.height)
            annotations.append(dict(hour=int(stage),category=category,stays=count,cohort_stays=n,
                panel_percent=float(row.panel_percent),percentage_label=row.percentage_label,
                display_text=label,placement=placement,within_node=bool(inside)))
            cursor+=count+gap
    # Stagger tiny adjacent-segment labels without altering node heights or counts.
    for stage_x in range(len(stages)):
        labels=[v for v in outside_labels if v[0]==stage_x]
        labels.sort(key=lambda v:ax.transData.transform(v[1].get_position())[1])
        previous=None
        for _,text,anchor_y in labels:
            px,py=ax.transData.transform(text.get_position());half=text.get_window_extent(renderer).height/2
            if previous is not None and py-half<previous+2:py=previous+2+half
            new_y=ax.transData.inverted().transform((px,py))[1]
            text.set_y(new_y);previous=py+half
            if abs(new_y-anchor_y)>1:
                tx=text.get_position()[0];edge=stage_x-bar/2 if tx<stage_x else stage_x+bar/2
                ax.plot([edge,tx],[anchor_y,new_y],color='#434A50',lw=.45,zorder=5)
    for x,(a,b) in enumerate(zip(stages[:-1],stages[1:])):
        rows=links[links.from_hour.eq(a)&links.to_hour.eq(b)].copy()
        from_order={key:i for i,key in enumerate(stage_categories[x])}
        to_order={key:i for i,key in enumerate(stage_categories[x+1])}
        rows['from_order']=rows.from_category.map(from_order);rows['to_order']=rows.to_category.map(to_order)
        rows=rows.sort_values(['from_order','to_order','origin'])
        source={key:0 for key in stage_categories[x]};target={key:0 for key in stage_categories[x+1]}
        for r in rows.itertuples(index=False):
            y0=positions[a,r.from_category][0]+source[r.from_category]
            y1=positions[b,r.to_category][0]+target[r.to_category]
            x0=x+bar/2;x1=x+1-bar/2;curve=.42*(x1-x0);h=int(r.stays)
            vertices=[(x0,y0),(x0+curve,y0),(x1-curve,y1),(x1,y1),(x1,y1+h),(x1-curve,y1+h),(x0+curve,y0+h),(x0,y0+h),(x0,y0)]
            codes=[MplPath.MOVETO,*([MplPath.CURVE4]*3),MplPath.LINETO,*([MplPath.CURVE4]*3),MplPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MplPath(vertices,codes),facecolor=styles[r.origin]['color'],alpha=.38,edgecolor='white',lw=.12,zorder=2))
            source[r.from_category]+=h;target[r.to_category]+=h
    check.update(dataset=str(nodes.dataset.iloc[0]),flow=str(nodes.flow.iloc[0]),horizon=str(nodes.horizon.iloc[0]),
        percentage_denominator='panel cohort_stays',percentage_decimals=1,node_percentage_labels=annotations,
        bar_width=bar,bar_width_mm=bar*ax.get_window_extent(renderer).width/(len(stages)-.28)*25.4/ax.figure.dpi)
    return check

def flow_styles():
    styles={f'Class {c}':dict(color=color,label=f'Class {c}') for c,color in COLORS.items()}
    styles.update({key:dict(value) for key,value in ENDPOINT_STYLES.items()})
    styles.update(Alive=dict(color='#173D56',label='Alive'),Died=dict(color='#71323C',label='Died'))
    return styles

def class_disposition():
    nodes=read('supp_disposition_nodes');links=read('supp_disposition_links');summary=read('supp_disposition_summary')
    styles=flow_styles();categories=[f'Class {c}' for c in COLORS]+['live_exit','icu_death','unknown']
    fig=plt.figure(figsize=(180/25.4,208/25.4));audits=[]
    for i,db in enumerate(DB):
        bottom=.64 if i==0 else .20
        ax=fig.add_axes([.10,bottom,.80,.25]);nn=nodes[nodes.dataset.eq(db)];ll=links[links.dataset.eq(db)]
        audits.append(exact_alluvial(ax,nn,ll,[1,3,5,7],[categories]*4,styles,label_mode='counts'))
        n=int(nn.cohort_stays.iloc[0])
        fig.text(.035,bottom+.29,chr(65+i),fontsize=10.5,fontweight='bold')
        fig.text(.08,bottom+.29,f'{DB[db]} · N={n:,}',fontsize=9.5,fontweight='bold')
        for stage,day in enumerate([1,3,5,7]):
            ax.text(stage,-.045,f'Day {day}',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=9,fontweight='bold')
            r=summary[summary.dataset.eq(db)&summary.day.eq(day)].iloc[0]
            line=f'ICU {int(r.resident):,}\nExit {int(r.live_icu_exit):,}\nDeath {int(r.icu_death):,}\nUnknown {int(r.unknown):,}'
            ax.text(stage,-.15,line,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8.5,linespacing=1.25)
    handles=[Patch(facecolor=styles[key]['color'],edgecolor=styles[key].get('edgecolor',styles[key]['color']),hatch=styles[key].get('hatch',''),label=styles[key].get('label',key)) for key in categories]
    fig.text(.5,.075,'Node labels: % of database N',ha='center',fontsize=8.5,color=GRAY)
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.012),ncol=3,fontsize=8.5,frameon=False)
    save(fig,'supp_figure_s6_class_disposition')
    records[-1]['flow_audit']=dict(type='fixed-class ICU disposition',checks=audits,
        days=[1,3,5,7],source='figure_source_data/supp_disposition_nodes.csv',
        link_source='figure_source_data/supp_disposition_links.csv',scaling='Each database independently scaled to its assigned N',
        class_switches=0,day7_resident=[1415,737],unknown_exit=[12,2])

def clinical_flows():
    dn,dl=read('supp_disposition_nodes'),read('supp_disposition_links')
    tn,tl=read('supp_treatment_flow_nodes'),read('supp_treatment_flow_links')
    definition=json.loads((S/'supplementary_flow_definitions.json').read_text())['definitions']['treatment_intensity']
    levels=definition['levels'];classes=[f'Class {c}' for c in COLORS];styles=flow_styles()
    for level,color in zip(levels,['#E9EDF1','#B6C0C9','#788A9B','#3A5067','#E6E8EB']):
        styles[level]=dict(color=color,label=level,hatch='////' if level=='Intensity unknown' else '')
    fig=plt.figure(figsize=(180/25.4,260/25.4));audits=[]
    for row,h in enumerate(['disposition','28','ICU']):
        for col,db in enumerate(DB):
            left=.055+.50*col;bottom=[.725,.435,.145][row]
            ax=fig.add_axes([left,bottom,.40,.195])
            if row==0:
                nn=dn[dn.dataset.eq(db)];ll=dl[dl.dataset.eq(db)];stages=[1,3,5,7]
                cats=classes+['live_exit','icu_death','unknown'];groups=[cats]*4;labels=['Day 1','Day 3','Day 5','Day 7']
            else:
                nn=tn[tn.dataset.eq(db)&tn.horizon.astype(str).eq(h)];ll=tl[tl.dataset.eq(db)&tl.horizon.astype(str).eq(h)];stages=[0,1,2]
                groups=[classes,levels,['Alive','Died','unknown'] if h=='28' else ['live_exit','icu_death','unknown']];labels=['Class','Treatment\nintensity','28-day\nmortality' if h=='28' else 'ICU\nmortality']
            assert int(nn.cohort_stays.iloc[0])==({'mimiciv':2698,'mimiciii':1327}[db])
            audits.append(exact_alluvial(ax,nn,ll,stages,groups,styles,label_mode='counts',bar_width=.50 if row==0 else .40,small_label_threshold=10))
            assert all(v['within_node'] for v in audits[-1]['node_percentage_labels'] if v['panel_percent']>=10),audits[-1]
            fig.text(left-.025,bottom+.225,chr(65+row*2+col),fontsize=10.5,fontweight='bold')
            fig.text(left+.015,bottom+.225,f'{DB[db]} · N={int(nn.cohort_stays.iloc[0]):,}',fontsize=9,fontweight='bold')
            for x,label in enumerate(labels):ax.text(x,-.05,label,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8.5,fontweight='bold')
    keys=[*classes,'unknown','Alive','Died','live_exit','icu_death']
    fig.legend(handles=[Patch(facecolor=styles[k]['color'],hatch=styles[k].get('hatch',''),label=styles[k].get('label',k)) for k in keys],loc='lower center',bbox_to_anchor=(.5,.066),ncol=4,frameon=False,fontsize=8.5)
    fig.legend(handles=[Patch(facecolor=styles[k]['color'],hatch=styles[k].get('hatch',''),label=k) for k in levels],loc='lower center',bbox_to_anchor=(.5,.018),ncol=3,frameon=False,fontsize=8.5)
    save(fig,'figure_3_class_treatment_outcomes')
    records[-1]['flow_audit']=dict(checks=audits,datasets=list(DB),horizons=['28','ICU'],panel_order=['IV disposition','III disposition','IV Day28','III Day28','IV ICU','III ICU'],ribbon_color='Persistent source trajectory class',definitions='figure_source_data/supplementary_flow_definitions.json')

def treatment_flows():
    clinical_flows()

def main():
    global OUTPUT, TREATMENT_SOURCE
    choices={'flow':flow,'management':management,'trajectories':trajectories,'phenotypes':phenotypes,'forest':mortality_forest,'survival':km,'subgroups':subgroup,'prediction':classifier,'supporting':supporting,'management-profiles':management_profiles,'treatment':treatment_forest,'treatment-probability':treatment_probability,'clinical-flows':clinical_flows,'treatment-flows':clinical_flows,'disposition':clinical_flows}
    parser = argparse.ArgumentParser()
    parser.add_argument('--figure', choices=['all', *choices], default='all')
    parser.add_argument('--treatment-source', choices=['released', 'refit'], default='released',
                        help='Use released S6/S7 aggregates or local outputs from fit_associations.R')
    args = parser.parse_args()
    TREATMENT_SOURCE = W / 'treatment/aggregate' if args.treatment_source == 'refit' else S
    with tempfile.TemporaryDirectory(prefix='critical-care-figures-') as temp:
        OUTPUT=Path(temp)
        if args.figure=='all':
            for function in dict.fromkeys(choices.values()):function()
        else:choices[args.figure]()
        for path in OUTPUT.iterdir():
            destination=P/'results/figures';destination.mkdir(parents=True,exist_ok=True)
            replacement=destination/(path.name+'.tmp');shutil.copy2(path,replacement);replacement.replace(destination/path.name)
        for record in records:
            if 'files' in record:record['files']=[str(HERE/Path(v).name) for v in record['files']]
        manifest=P/'results/figure_export_manifest.json'
        if args.figure!='all' and manifest.exists():
            previous=json.loads(manifest.read_text());changed={r['stem'] for r in records}
            obsolete={'figure_3_early_management_profiles','supp_figure_s6_class_disposition','supp_figure_s7_treatment_intensity_outcomes'}
            records[:]=[r for r in previous if r['stem'] not in changed|obsolete]+records
        manifest.write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
    print(f'Exported {len(records)} figure groups as PDF, native 300 dpi RGB TIFF and document PNG.')

if __name__=='__main__':main()
