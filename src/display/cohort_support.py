"""Local cohort arithmetic and measured text wrapping for the study renderer."""
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.font_manager import FontProperties
PALETTES={'primary':{'header':'#F7CFAD','fill':'#FBE8D8'},'external':{'header':'#B8C4D3','fill':'#E2E7EE'}}
def node_text(node):
    return f"{node['label']}\n(n={node['n']:,})"
def validate_config(config):
    steps=[]
    for cohort in config['cohorts']:
        for previous,current,excluded in zip(cohort['nodes'],cohort['nodes'][1:],cohort['exclusions']):
            if previous['n']-excluded['n']!=current['n']:raise ValueError('Cohort counts do not reconcile')
            steps.append({'excluded_n':excluded['n']})
    return config,steps
def wrap_text(text,width_mm,renderer,font,size,bold):
    prop=FontProperties(family=font,size=size,weight='bold' if bold else 'normal')
    pixels=width_mm/25.4*renderer.dpi
    lines=[]
    for paragraph in text.split('\n'):
        line=''
        for word in paragraph.split():
            candidate=(line+' '+word).strip()
            if line and renderer.get_text_width_height_descent(candidate,prop,False)[0]>pixels:
                lines.append(line);line=word
            else:line=candidate
        lines.append(line)
    return '\n'.join(lines)
