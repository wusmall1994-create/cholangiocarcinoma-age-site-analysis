"""Redraw main figures from existing aggregate estimates only; no model fitting.

Run: python cancer_medicine_figures.py --results RESULTS --out OUTPUT
Lines connect reported ages only, not newly predicted ages. Exports are 600-dpi TIFF.
"""
from pathlib import Path
import argparse, csv
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

parser=argparse.ArgumentParser()
parser.add_argument('--results',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
args.out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.7})
BLUE='#167b8c';GRAY='#505b66';ORANGE='#b66b23'
def rows(name):
    with (args.results/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def save(fig,n):
    fig.savefig(args.out/f'Figure_{n}.tiff',dpi=600,format='tiff',pil_kwargs={'compression':'tiff_lzw'},facecolor='white')
    plt.close(fig)

# Figure 1: distinguish all-age counts from exact-age models.
fig,ax=plt.subplots(figsize=(7.2,6.3));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
def box(x,y,w,h,text,color='#edf3f7'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008',facecolor=color,edgecolor='#48606e',linewidth=.8))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=9,linespacing=1.4)
box(.10,.85,.8,.12,'SEER 17 malignant bile-duct records, 2000–2023: 45,243\nHistology 8160/8162, ages 15 or older: 31,086')
box(.10,.67,.8,.12,'Diagnosed 2004–2023: 28,654 records / 28,605 patients\nPrimary exact-age, known-cause analysis: 25,400 records')
box(.02,.42,.45,.16,'Original binary risk analysis\n2004–2018: 15,578 records\nUnpenalized cause-specific models\nOwn cohort as standardization target')
box(.53,.42,.45,.16,'Harmonized definition analysis\n2010–2018: 11,400 records\nFour binary definitions (D1–D4)\nFixed weak ridge in both hazards')
box(.53,.16,.45,.18,'Three-category EOD analysis\n11,361 records / 11,347 patients\n8160 fitting cohort: 10,882 records\nCommon EOD target for all six analyses','#e6f1ee')
for start,end in [((.5,.85),(.5,.79)),((.32,.67),(.25,.58)),((.68,.67),(.75,.58)),((.755,.42),(.755,.34))]:
    ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'->','color':'#48606e'})
ax.text(.02,.28,'Main exact-age models: ages 15–84;\ncancer death with other death competing.\nModel and standardization targets\nare reported separately.',fontsize=9,va='top')
ax.text(.5,.055,'Numbers denote tumor records unless patients are explicitly stated.\nThe harmonized extension shares 1,000 patient resamples across definitions.',ha='center',fontsize=8)
fig.tight_layout();save(fig,1)

# Figure 3: no causal or equivalence interpretation of components or intervals.
base=Path('composition_decomposition/paired_1000')
est=rows(base/'estimates_with_ci.csv');paired=rows(base/'paired_differences_with_ci.csv')
def get(data,a,c):return next(r for r in data if r['analysis']==a and int(r['horizon_months'])==60 and r['component']==c)
def point(ax,r,y,color):
    x=float(r['estimate'])*100;lo=float(r['lower95'])*100;hi=float(r['upper95'])*100
    ax.errorbar(x,y,xerr=[[x-lo],[hi-x]],fmt='o',color=color,elinewidth=1.1,capsize=3,markersize=4)
fig,axes=plt.subplots(3,1,figsize=(7.2,9.1),gridspec_kw={'height_ratios':[1.25,1.15,1]},layout='constrained')
labels={'D1_topography':'D1  Topography','D2_exclude_discordant':'D2  Exclude C22.1 / 8162','D3_8160_only':'D3  Histology 8160 only','D4_reallocate_8162':'D4  Hypothetical reallocation','EOD':'EOD  Synthetic mixture','EOD_8160':'EOD  Histology 8160 only'}
ax=axes[0]
for i,a in enumerate(labels):point(ax,get(est,a,'total'),i,BLUE if a.startswith('EOD') else GRAY)
ax.set_yticks(range(6),labels.values());ax.invert_yaxis();ax.set_title('a  Age 75 minus age 50: change in five-year risk difference',loc='left',fontweight='bold');ax.set_xlabel('Change in eCCA minus iCCA risk difference (percentage points)')
components={'composition':'Anatomical composition','perihilar_risk':'Perihilar risk pattern','distal_risk':'Distal risk pattern','iCCA_reference':'iCCA reference risk','total':'Total change'}
ax=axes[1]
for i,c in enumerate(components):
    for a,off,color in [('EOD',-.13,BLUE),('EOD_8160',.13,ORANGE)]:point(ax,get(est,a,c),i+off,color)
ax.set_yticks(range(5),components.values());ax.invert_yaxis();ax.set_title('b  Descriptive decomposition of the synthetic mixture contrast',loc='left',fontweight='bold');ax.set_xlabel('Descriptive component (percentage points)')
ax.plot([],[],'o',color=BLUE,label='EOD');ax.plot([],[],'o',color=ORANGE,label='EOD: 8160 only');ax.legend(frameon=False,fontsize=8,loc='upper left')
ax=axes[2]
names=['D2_exclude_discordant_minus_D1_topography','D3_8160_only_minus_D1_topography','D4_reallocate_8162_minus_D1_topography','EOD_8160_minus_EOD']
for i,a in enumerate(names):point(ax,get(paired,a,'total'),i,GRAY)
ax.set_yticks(range(4),['D2 − D1','D3 − D1','D4 − D1','EOD 8160 − EOD']);ax.invert_yaxis();ax.set_title('c  Paired differences in the age-related change',loc='left',fontweight='bold');ax.set_xlabel('Difference in age-related change (percentage points)')
for ax in axes:ax.axvline(0,color='#777777',ls=':',lw=.8);ax.grid(axis='x',alpha=.15)
save(fig,3)
print('Figures 1 and 3 generated from aggregate estimates.')
