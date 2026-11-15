"""Study-local matplotlib export and aggregate count validation."""
from pathlib import Path
import numpy as np
import matplotlib as mpl
ENDPOINT_STYLES={
'pre_icu':dict(color='#FFFFFF',edgecolor='#68717A',hatch='..',label='Not yet in ICU'),
'live_exit':dict(color='#173D56',edgecolor='#173D56',hatch='',label='Live ICU exit'),
'icu_death':dict(color='#71323C',edgecolor='#71323C',hatch='',label='ICU death'),
'unknown':dict(color='#E6E8EB',edgecolor='#747B83',hatch='////',label='Unknown state / end')}
def apply_style():
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':9,'axes.labelsize':8.5,'axes.titlesize':9,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,'xtick.direction':'out','ytick.direction':'out','xtick.major.width':.6,'ytick.major.width':.6,'lines.linewidth':1,'hatch.linewidth':.4,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.dpi':150,'savefig.dpi':300,'savefig.bbox':None,'axes.unicode_minus':False})
def export_figure(fig,output_stem):
    stem=Path(output_stem);files=[]
    for ext in ['pdf','svg','tif']:
        path=stem.with_suffix('.'+ext)
        options={'dpi':300,'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tif' else {}
        fig.savefig(path,facecolor='white',bbox_inches=None,**options);files.append(str(path))
    return {'files':files,'width_mm':float(fig.get_figwidth()*25.4),'height_mm':float(fig.get_figheight()*25.4),'dpi':300}
def validate_flow(nodes,links,hours,categories=None):
    for hour in hours:
        if not np.isclose(nodes.loc[nodes.hour.eq(hour),'stays'].sum(),nodes.loc[nodes.hour.eq(hours[0]),'stays'].sum()):raise ValueError('Flow population is not conserved')
    for left,right in zip(hours[:-1],hours[1:]):
        segment=links[links.from_hour.eq(left)&links.to_hour.eq(right)]
        for hour,column in [(left,'from_category'),(right,'to_category')]:
            expected=nodes[nodes.hour.eq(hour)].groupby('category').stays.sum()
            actual=segment.groupby(column).stays.sum().reindex(expected.index,fill_value=0)
            if not np.array_equal(expected.to_numpy(),actual.to_numpy()):raise ValueError('Flow edge count mismatch')
    return {'conserved':True,'cohort_stays':int(nodes.loc[nodes.hour.eq(hours[0]),'stays'].sum())}
