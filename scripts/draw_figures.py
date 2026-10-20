import os

from pathlib import Path

import csv,json

import matplotlib as mpl

import matplotlib.pyplot as plt

from matplotlib.patches import FancyBboxPatch



O=Path(os.environ['CCA_FIGURE_OUTPUT_DIR']);Q=Path(os.environ['CCA_RESULTS_DIR']);O.mkdir(parents=True,exist_ok=True)

mpl.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.6,'xtick.labelsize':8,'ytick.labelsize':8,'savefig.facecolor':'white'})

contract={'backend':'Python matplotlib; matplotlib','format':'TIFF only, 600 dpi, LZW; internal previews TIFF','Fig2':{'claim':'A recorded death can be assigned differently by endpoint rules while remaining an all-cause event.','archetype':'classification schematic','evidence':'Official rules plus observed 621-record crosswalk; arrows describe classification only','risk':'No causal diagram, error adjudication, or diagnostic gold standard'},'Fig3':{'claim':'Age-specific site contrasts are sensitive to endpoints; sequence strata have distinct targets.','archetype':'quantitative grid','uncertainty':'Pointwise percentile 95% CIs, fixed theta=.1, 2000 paired patient-cluster samples in unified G','risk':'Both common endpoints shown; different x scales explicitly labeled; no between-stratum causal decomposition'}}

(Q/'figure_contract.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')



# Cohort map prioritizes the endpoint and sequence analysis populations.

f=plt.figure(figsize=(7.2,5.3));aa=f.add_axes([0,0,1,1]);aa.set_axis_off();aa.set_xlim(0,1);aa.set_ylim(0,1)

def cb(x,y,w,h,t,size=8.5,fill='#F2F5F6'):

 aa.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.008,rounding_size=.007',lw=.8,edgecolor='#687985',facecolor=fill))

 aa.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=size,linespacing=1.3)

def ca(a,b):aa.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','lw':.8,'color':'#7c858c','shrinkA':2,'shrinkB':2})

cb(.17,.88,.66,.085,'SEER 17 malignant C22.1 / C24.0 records, 2000–2023\n45,243 records',9.3)

ca((.5,.875),(.5,.84))

cb(.17,.725,.66,.11,'Histology 8160 / 8162, age 15 years or older\n31,086 records in 2000–2023\n28,654 records in 2004–2023',9)

ca((.4,.72),(.17,.655));ca((.61,.72),(.665,.655))

cb(.025,.47,.29,.18,'Original hazard hypothesis\n2004–2023, exact age, known cause\n25,400 records\n25,353 patients',8.5)

cb(.355,.49,.62,.16,'Five-year risk population, 2004–2018, ages 15–84\n15,780 records / 15,763 patients\n202 unknown-cause deaths included in direct all-cause models',8.8)

ca((.665,.485),(.665,.435))

cb(.355,.285,.62,.145,'Known-cause endpoint and sequence analyses\n15,578 records / 15,561 patients\nEndpoints share the same within-population target and resamples',8.8,fill='#ECF2F3')

ca((.54,.28),(.50,.23));ca((.80,.28),(.835,.23))

cb(.355,.105,.29,.12,'Registry-recorded single primary\n12,222 records and patients',8.3)

cb(.685,.105,.29,.12,'Multiple-primary sequence\n3,356 records / 3,339 patients',8.3)

aa.text(.025,.345,'Separate anatomical extension',fontsize=8.5,fontweight='bold')

aa.text(.025,.30,'2010–2018 EOD target',fontsize=8.3)

aa.text(.025,.26,'11,361 classifiable records',fontsize=8.3)

aa.text(.025,.22,'Original SEER death endpoint',fontsize=8.3)

aa.text(.025,.18,'Detailed eligibility in Table S18',fontsize=8.3)

aa.text(.665,.035,'Sequence strata use their own empirical targets',ha='center',fontsize=8,color='#555f66')

f.savefig(O/'Fig1.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(f)



fig=plt.figure(figsize=(7.2,4.9));ax=fig.add_axes([0,0,1,1]);ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)

colors=['#245B78','#AB742A','#68747D']

def box(x,y,w,h,text,color='#49555c',fill='white',size=8.5):

    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008,rounding_size=0.008',lw=.8,edgecolor=color,facecolor=fill))

    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=size,color='#20282d',linespacing=1.35)

def arrow(a,b,color='#7c858c'):

    ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':color,'lw':.85,'shrinkA':2,'shrinkB':2})

box(.31,.84,.38,.105,'Same recorded death\nUnderlying COD label',fill='#F2F5F6',size=10)

xs=[.025,.355,.685];w=.29

for x in xs:arrow((.5,.83),(x+w/2,.805))

texts=[('SEER attribution','Primary site + primary sequence\n+ underlying cause of death','Attributed endpoint\nor complementary death'),('Common biliary set','Intrahepatic Bile Duct\nor Other Biliary','Set member: endpoint\nOther known death: complement'),('Common hepatobiliary set','Biliary set + Liver\n+ Gallbladder','Set member: endpoint\nOther known death: complement')]

for x,col,(head,rule,outcome) in zip(xs,colors,texts):

    ax.text(x+w/2,.78,head,ha='center',va='bottom',fontweight='bold',fontsize=9.2,color=col)

    box(x,.58,w,.17,rule,color=col)

    arrow((x+w/2,.575),(x+w/2,.52),col)

    box(x,.405,w,.11,outcome,color=col,size=8.2)

ax.text(.03,.345,'Observed example',fontweight='bold',fontsize=9.5)

ax.text(.03,.298,'621 eCCA records in multiple-primary sequences, COD label Intrahepatic Bile Duct',fontsize=8.6)

for x,col,txt in zip(xs,colors,['Complementary death','Endpoint death','Endpoint death']):box(x,.205,w,.065,txt,color=col,fill='#F7F8F8',size=9)

box(.21,.085,.58,.065,'Direct all-cause endpoint: all 621 are death events',fill='#F2F5F6',size=9)

ax.text(.5,.025,'Classification paths, not causal mechanisms or adjudicated errors',ha='center',fontsize=8,color='#555f66')

fig.savefig(O/'Fig2.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)



with (Q/'ridge_contrasts_with_CI.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))

defs=['SEER_original','Common_biliary','Common_hepatobiliary'];labs=['SEER attribution','Common biliary','Common hepatobiliary']

groups=[('All records',15578,'All'),('Registry-recorded\nsingle primary',12222,'Single_primary'),('Multiple-primary\nsequence',3356,'Multiple_primary_sequence')]

fig,axs=plt.subplots(2,3,figsize=(7.2,5.5),sharey=True);source=[]

for col,(lab,n,g) in enumerate(groups):

 for row,c in enumerate(['RD_age75','RD_change_75_minus_50']):

  ax=axs[row,col]

  for j,(de,color) in enumerate(zip(defs,colors)):

   x=next(x for x in rows if x['sequence_group']==g and x['definition']==de and x['state']=='cancer' and x['contrast']==c)

   v,lo,hi=[100*float(x[k]) for k in ['estimate','lower95','upper95']]

   ax.errorbar(v,2-j,xerr=[[v-lo],[hi-v]],fmt='o',color=color,markersize=4,capsize=2,lw=1)

   if col==1:
    label=f'{v:.2f} ({lo:.2f} to {hi:.2f})'.replace('-', '−')
    ax.text(.03,2-j-.29,label,transform=ax.get_yaxis_transform(),fontsize=6.8,ha='left',va='top',color=color)
   source.append(dict(group=g,definition=de,contrast=c,estimate_pp=v,lower95_pp=lo,upper95_pp=hi))

  ax.axvline(0,color='#999999',lw=.7,ls='--',zorder=0)

  ax.set_yticks([2,1,0],labs);ax.set_ylim(-.55,2.65)

  ax.set_xlim((-65,20) if row==0 else (-25,20));ax.set_xticks([-60,-40,-20,0,20] if row==0 else [-20,-10,0,10,20])

  ax.set_xlabel('RD at age 75 (pp)' if row==0 else 'RD75 − RD50 (pp)')

  ax.text(-.09,1.05,'abcdef'[row*3+col],transform=ax.transAxes,fontweight='bold',fontsize=10)

  if row==0:ax.set_title(f'{lab}\nN = {n:,}',fontsize=8.5,pad=10)

# Paired contrast uses its own bootstrap interval, not marginal interval subtraction.
pair=next(x for x in rows if x['sequence_group']=='All' and x['definition']=='Common_biliary_minus_SEER' and x['state']=='cancer' and x['contrast']=='RD_change_75_minus_50')
pv,pl,ph=[100*float(pair[k]) for k in ['estimate','lower95','upper95']]
a=axs[1,0]
a.annotate('',xy=(16,1),xytext=(16,2),arrowprops={'arrowstyle':'|-|','lw':.8,'color':colors[1]})
a.text(13.5,1.5,f'{pv:.2f} pp',rotation=90,ha='center',va='center',fontsize=7,color=colors[1])
fig.text(.235,.022,f'Paired common biliary − SEER age contrast: {pv:.2f} pp (95% CI {pl:.2f} to {ph:.2f})',fontsize=7.5,color=colors[1])
fig.subplots_adjust(left=.235,right=.98,top=.82,bottom=.19,wspace=.18,hspace=.65)

fig.savefig(O/'Fig3.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)

with (Q/'Fig3_source.csv').open('w',encoding='utf-8',newline='') as f:

 w=csv.DictWriter(f,fieldnames=source[0]);w.writeheader();w.writerows(source)

print('Saved three TIFF figures')
