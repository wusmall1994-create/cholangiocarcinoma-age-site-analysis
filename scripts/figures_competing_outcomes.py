# Supplementary analysis plots. For manuscript Figures 1–3 use draw_figures.py.
"""Quantitative grids: competing-outcome distinction and exploratory model sensitivity.
Uses aggregate estimates only. TIFF output, 600 dpi; pointwise intervals.
Usage: python figures_competing_outcomes.py --clinical results/clinical_extension --results results --out figures
"""
from pathlib import Path
import argparse,csv,shutil
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--clinical',type=Path,default=Path('results/clinical_extension'))
parser.add_argument('--results',type=Path,default=Path('results'))
parser.add_argument('--out',type=Path,default=Path('figures'))
a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=True)

def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def pick(data,**kw):
    z=[r for r in data if all(str(r[k])==str(v) for k,v in kw.items())];assert len(z)==1,(kw,len(z));return z[0]
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'legend.frameon':False})
BLUE='#167b8c';GRAY='#505b66';ORANGE='#b66b23';PURPLE='#705e9b'
def save(fig,name):
    fig.savefig(a.out/name,dpi=600,format='tiff',pil_kwargs={'compression':'tiff_lzw'},facecolor='white')
    plt.close(fig)
def err(ax,z,x,color,marker,label,scale=100):
    y=np.array([float(r['estimate'])*scale for r in z]);lo=np.array([float(r['lower95'])*scale for r in z]);hi=np.array([float(r['upper95'])*scale for r in z])
    ax.errorbar(x,y,yerr=[y-lo,hi-y],fmt=marker+'-',color=color,capsize=3,lw=1.3,markersize=4.5,label=label)

st=rows(a.clinical/'original_states_with_CI.csv');rd=rows(a.clinical/'original_RDs_with_CI.csv')
ages=[40,50,65,75]
fig,axes=plt.subplots(2,1,figsize=(7.2,7.0),layout='constrained',gridspec_kw={'height_ratios':[1,1.25]})
for site,col,mark in [('iCCA',GRAY,'s'),('eCCA',BLUE,'o')]:
    z=[pick(st,horizon_months=60,age=age,exposure_level=site,state='cancer') for age in ages]
    err(axes[0],z,ages,col,mark,site)
axes[0].set_title('a  Five-year standardized cancer-death risks',loc='left',fontweight='bold')
axes[0].set(ylabel='Cancer-death risk (%)',ylim=(60,89));axes[0].legend(loc='lower left',ncol=2)
for state,col,mark,label in [('cancer',BLUE,'o','Cancer death'),('other',ORANGE,'^','Other death'),('all_death',PURPLE,'s','All-cause death')]:
    z=[pick(rd,horizon_months=60,age=age,state=state) for age in ages]
    err(axes[1],z,ages,col,mark,label)
axes[1].set_title('b  Competing outcomes and all-cause death',loc='left',fontweight='bold')
axes[1].set(ylabel='eCCA minus iCCA (percentage points)',xlabel='Age at diagnosis, years',ylim=(-21,21))
axes[1].axhline(0,color='#888888',lw=.7);axes[1].legend(loc='lower left',ncol=3,fontsize=8)
axes[1].text(.02,.91,'At age 75: all-cause RD +1.08 pp (95% CI −0.09 to 2.22)',transform=axes[1].transAxes,fontsize=8)
for ax in axes:ax.set_xticks(ages);ax.set_xlim(37,78);ax.grid(axis='y',alpha=.15)
save(fig,'Figure_2.tiff')

# S3 displays composition and relative effects.
counts=rows(a.results/'tables/diagnostic_schema3_counts.csv')
binary=rows(a.results/'tables/diagnostic_exact_age_effects.csv')
sub=rows(a.results/'tables/diagnostic_schema3_effect_curves.csv')
fig,axes=plt.subplots(3,1,figsize=(7.2,7.5),layout='constrained')
bands=['15-39','40-64','65-84'];frac=[];ns=[]
for band in bands:
    z=[r for r in counts if r['age_group_exact']==band and r['schema3']!='iCCA'];den=sum(int(r['n']) for r in z)
    frac.append(sum(int(r['n']) for r in z if r['schema3']=='Distal eCCA')/den);ns.append(den)
x=np.arange(3);f=np.array(frac)*100
axes[0].bar(x,100-f,color=GRAY,label='Perihilar');axes[0].bar(x,f,bottom=100-f,color=BLUE,label='Distal')
for xx,v,n in zip(x,f,ns):axes[0].text(xx,100-v/2,f'{v:.1f}%',ha='center',va='center',color='white',fontsize=9)
axes[0].set_xticks(x,[f'{b.replace("-","–")}\nN={n:,}' for b,n in zip(bands,ns)])
axes[0].set(ylabel='eCCA composition (%)',ylim=(0,118));axes[0].set_yticks([0,20,40,60,80,100]);axes[0].legend(ncol=2,loc='upper right',fontsize=8)
axes[0].set_title('a  Detailed anatomy by age group, 2010–2023',loc='left',fontweight='bold')
for ax,data,contrasts,title in [(axes[1],binary,['eCCA'],'b  Binary Model 1, 2004–2023'),(axes[2],sub,['Perihilar eCCA','Distal eCCA'],'c  Detailed anatomy Model 1, 2010–2023')]:
    for site,col,marker in zip(contrasts,[GRAY,BLUE],['s','o']):
        z=[pick(data,age=ag,model='Model 1',contrast=site) for ag in ages]
        err(ax,z,ages,col,marker,site,scale=1)
    ax.axhline(1,color='#888888',ls=':',lw=.8);ax.set(ylabel='csHR versus iCCA',xlabel='Age at diagnosis, years',ylim=(.45,1.25));ax.set_xticks(ages);ax.legend(ncol=2,loc='lower left');ax.set_title(title,loc='left',fontweight='bold')
save(fig,'Supplementary_Figure_S3.tiff')

dec=rows(a.clinical/'harmonized_decomposition_with_CI.csv');paired=rows(a.clinical/'interval_minus_PH_paired_CI.csv')
fig,axes=plt.subplots(2,1,figsize=(7.2,6.5),layout='constrained',gridspec_kw={'height_ratios':[1.6,1]})
comps={'composition':'Composition','perihilar_risk':'Perihilar risk','distal_risk':'Distal risk','iCCA_reference':'iCCA reference','total':'Total change'}
def horizontal(ax,r,y,col,marker):
    e=float(r['estimate'])*100;l=float(r['lower95'])*100;u=float(r['upper95'])*100
    ax.errorbar(e,y,xerr=[[e-l],[u-e]],fmt=marker,color=col,capsize=3,markersize=5)
for i,c in enumerate(comps):
    for model,offset,col,mark in [('PH',-.12,GRAY,'s'),('interval',.12,BLUE,'o')]:horizontal(axes[0],pick(dec,horizon_months=60,state='cancer',component=c,model=model),i+offset,col,mark)
axes[0].set_yticks(range(5),comps.values());axes[0].invert_yaxis();axes[0].set_title('a  Cancer-death RD change and descriptive components',loc='left',fontweight='bold')
axes[0].plot([],[],'s',color=GRAY,label='PH');axes[0].plot([],[],'o',color=BLUE,label='Interval-specific');axes[0].legend(loc='lower right')
for i,stt in enumerate(['cancer','other','all_death']):horizontal(axes[1],pick(paired,horizon_months=60,state=stt,component='total'),i,PURPLE,'o')
axes[1].set_yticks(range(3),['Cancer death','Other death','All-cause death']);axes[1].invert_yaxis();axes[1].set_title('b  Paired interval-minus-PH differences in total change',loc='left',fontweight='bold')
for ax in axes:ax.axvline(0,color='#888888',ls=':',lw=.7);ax.set_xlabel('Percentage points');ax.grid(axis='x',alpha=.15)
save(fig,'Supplementary_Figure_S4.tiff')
print('Figure 2 and Supplementary Figures S3 and S4 generated.')
