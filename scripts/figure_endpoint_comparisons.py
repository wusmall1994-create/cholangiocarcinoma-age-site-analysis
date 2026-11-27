import csv
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
here = Path(__file__).resolve().parents[1]
stage = here/'figures'
stage.mkdir(exist_ok=True)
source = here/'results'/'endpoint_sensitivity'
with (source/'key_results_percentage_points.csv').open(newline='') as f:
    data = {(r['sequence_group'],r['definition'],r['contrast']):
            tuple(float(r[k]) for k in ('estimate','lower95','upper95')) for r in csv.DictReader(f)}

plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.spines.left':False,'axes.linewidth':.7,
                     'pdf.fonttype':42})
blue, gray, green, red = '#17546b','#4d5c67','#447350','#954525'
def signed(x):
    return f'{x:+.2f}'.replace('-', '−')

def forest(ax, rows, ylim, xlim, step, title, xlabel, top):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(MultipleLocator(step))
    ax.axvline(0, color='#949aa0', linestyle='--', linewidth=.8, zorder=0)
    ax.grid(axis='x',color='#eceef0',linewidth=.6,zorder=0)
    ax.set_xlabel(xlabel)
    ax.set_title(title,loc='left',pad=12,fontweight='bold',fontsize=10)
    ax.text(1.065,top,'Estimate (95% CI)',transform=ax.get_yaxis_transform(),fontweight='bold',va='center')
    for y, label, key, color in rows:
        v,lo,hi = data[key]
        ax.errorbar(v,y,xerr=[[v-lo],[hi-v]],fmt='o',color=color,capsize=2.5,linewidth=1.2,markersize=4.6)
        ax.text(-.045,y,label,transform=ax.get_yaxis_transform(),va='center',ha='right',color='#242a2e')
        ax.text(1.065,y,f'{signed(v)} ({signed(lo)} to {signed(hi)})',transform=ax.get_yaxis_transform(),va='center',color=color)

fig=plt.figure(figsize=(8.6,7.4))
g=fig.add_gridspec(2,1,left=.39,right=.755,bottom=.09,top=.92,hspace=.60,height_ratios=[3.15,1.85])
a,b=(fig.add_subplot(g[i]) for i in range(2))
age='RD_age75'
forest(a,[
    (7,'SEER attribution',('All','SEER_original',age),blue),
    (6,'Common biliary',('All','Common_biliary',age),blue),
    (5,'Common hepatobiliary',('All','Common_hepatobiliary',age),blue),
    (4,'621-death reclassification',('All','Scenario621',age),blue),
    (2,'SEER attribution',('Single_primary','SEER_original',age),gray),
    (0,'Direct all-cause',('Expanded','Direct_allcause_including_unknown',age),green),
],(-.55,9),(-20,10),5,'A   Site difference at age 75','Five-year risk difference (percentage points)',8.45)
for y,label in [(8,'Cause-known cohort (N = 15,578)'),(3,'Single-primary subset (N = 12,222)'),(1,'Expanded cohort (N = 15,780)')]:
    a.text(-.045,y,label,transform=a.get_yaxis_transform(),ha='right',va='center',fontweight='bold',fontsize=8.5)
change='RD_change_75_minus_50'
with (here/'results/sequence_sensitivity/A40_summary.csv').open() as f:
    row = next(r for r in csv.DictReader(f) if r['definition']=='Common_biliary_minus_SEER' and r['contrast']=='RD_age75')
    data[('First','Common_biliary_minus_SEER',age)] = tuple(100*float(row[k]) for k in ('estimate','lower95','upper95'))
forest(b,[
    (4,'Full cohort, age 75',('All','Common_biliary_minus_SEER','RD_age75'),blue),
    (3,'First-sequenced subset, age 75',('First','Common_biliary_minus_SEER',age),green),
    (2,'Full cohort, age 50',('All','Common_biliary_minus_SEER','RD_age50'),blue),
    (1,'Common biliary vs SEER, 75 minus 50',('All','Common_biliary_minus_SEER',change),blue),
    (0,'621-death scenario vs SEER, 75 minus 50',('All','Scenario621_minus_SEER',change),gray),
],(-.55,5.3),(-2,21),5,'B   Paired effect of changing the endpoint','Change in five-year risk difference (percentage points)',5.0)
fig.canvas.draw()

def export(fig,n):
    fig.savefig(stage/f'Figure_{n}.tif',dpi=600,pil_kwargs={'compression':'tiff_lzw'})

export(fig,2)
plt.close(fig)

with (here/'results/sequence_sensitivity/Fig3_counts.csv').open() as f:
    crosswalk={(r['sequence'],r['site'],r['endpoint']):(int(r['seer_only']),int(r['common_only'])) for r in csv.DictReader(f)}
assert crosswalk['Multiple','eCCA','Biliary']==(73,590)

fig, ax=plt.subplots(figsize=(7.5,3.45))
fig.subplots_adjust(left=.06,right=.94,bottom=.06,top=.92)
ax.set_axis_off()
ax.set(xlim=(0,1),ylim=(0,1))
ax.text(0,.96,'Deaths with different endpoint classifications',fontsize=12,fontweight='bold')
ax.text(0,.865,'Deaths within 60 months of diagnosis',color=gray)
ax.text(.45,.74,'SEER vs common biliary',ha='center',fontweight='bold')
ax.text(.825,.74,'SEER vs common hepatobiliary',ha='center',fontweight='bold')
ax.text(.02,.63,'Primary sequence and site',fontweight='bold')
positions=(.365,.535,.74,.91)
for x,label in zip(positions,['SEER only','Common only','SEER only','Common only']):
    ax.text(x,.63,label,ha='center',fontweight='bold',fontsize=8)
for (seq,site),y in zip([('Single','iCCA'),('Single','eCCA'),('Multiple','iCCA'),('Multiple','eCCA')],[.51,.40,.29,.18]):
    counts=crosswalk[seq,site,'Biliary']+crosswalk[seq,site,'Hepatobiliary']
    ax.hlines(y-.06,.02,.98,color='#e1e6e8',lw=.65)
    ax.text(.02,y,f'{seq}-primary, {site}',va='center')
    for i,(x,count) in enumerate(zip(positions,counts)):
        focus=i==1 and count==590
        ax.text(x,y,f'{count:,}',ha='center',va='center',color=red if focus else (blue if i%2 else gray),fontweight='bold' if focus else 'normal')
ax.vlines(.637,.12,.80,color='#e1e6e8',lw=.8)
ax.text(.02,.04,'590 of the 621 targeted deaths identified across full follow-up occurred within 60 months.',
        color=red,fontsize=8.5,va='center')
fig.canvas.draw()
export(fig,3)

# Main analysis populations. Counts are tumor records.
fig, ax = plt.subplots(figsize=(8,5))
ax.set(xlim=(0,1),ylim=(0,1))
ax.axis('off')
def box(x,y,label):
    ax.text(x,y,label,ha='center',va='center',fontsize=10,
            bbox=dict(boxstyle='round,pad=0.65',facecolor='#f2f7f8',edgecolor='#17546b'))
def arrow(x,y,u,v):
    ax.annotate('',xy=(u,v),xytext=(x,y),arrowprops=dict(arrowstyle='->',color=gray))
box(.5,.90,'Malignant C22.1 or C24.0 tumors, 2000–2023\n45,243 records after export exclusions')
box(.5,.65,'Eligible histology: 31,086 records\nDiagnosed in 2004–2023: 28,654')
arrow(.5,.82,.5,.74)
box(.24,.39,'Age-by-site hazard analysis\n2004–2023, ages 15–84\n25,400 cause-known records')
box(.76,.39,'Direct all-cause mortality\n2004–2018, ages 15–84\n15,780 records')
arrow(.45,.57,.25,.50)
arrow(.55,.57,.75,.50)
box(.76,.12,'Paired endpoint comparison\n15,578 cause-known records')
arrow(.76,.28,.76,.20)
ax.text(.70,.235,'Exclude 202 deaths of unknown cause',ha='right',va='center',fontsize=9)
export(fig,1)
plt.close(fig)
