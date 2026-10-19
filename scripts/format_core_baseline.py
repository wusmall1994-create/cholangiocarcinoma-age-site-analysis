from pathlib import Path
import csv,json,os
Q=Path(os.environ['CCA_RESULTS_DIR'])
def read(n):return list(csv.DictReader((Q/n).open(encoding='utf-8-sig')))
a=read('age_summary.csv');ns={x['site_group']:int(x['N']) for x in a};counts=read('baseline_counts.csv')
def n(v,l,s):return sum(int(x['n']) for x in counts if x['variable']==v and x['level'] in l and x['site']==s)
def fmt(k,s):return f'{k:,} ({100*k/ns[s]:.1f}%)' if k==0 or 100*k/ns[s]>=.05 else f'{k:,} (<0.1%)'
rows=[['Characteristic','Level',f'iCCA (N={ns["iCCA"]:,})',f'eCCA (N={ns["eCCA"]:,})'],['Age, years','Median (IQR)']+[next(f'{float(x["median"]):g} ({float(x["q1"]):g}–{float(x["q3"]):g})' for x in a if x['site_group']==s) for s in ['iCCA','eCCA']]]
for i,l in enumerate(['15–39','40–49','50–64','65–74','75–84']):rows.append(['Age group' if i==0 else '',l]+[fmt(sum(int(x['n']) for x in read('age_counts.csv') if x['site']==s and x['level']==l),s) for s in ['iCCA','eCCA']])
spec=[('Sex','sex_factor',['Female','Male']),('Race/ethnicity','race_factor',['Non-Hispanic White','Non-Hispanic Black','Non-Hispanic Asian or Pacific Islander','Non-Hispanic American Indian/Alaska Native','Hispanic (All Races)','Non-Hispanic Unknown Race']),('Marital status','marital4',['Married','Unmarried','Widowed','Unknown']),('County household income, US$','income3',['<70,000','70,000-99,999','>=100,000','Unknown']),('County rurality','rural3',['Metropolitan','Nonmetropolitan','Unknown']),('Diagnosis era','era4',['2004-2009','2010-2015','2016-2019']),('Summary stage','stage4',['Localized','Regional','Distant','Unknown']),('Grade','grade4',['I','II','III','Unknown']),('Primary-site surgery','surgery_primary',['Yes'])]
for label,v,ls in spec:
 for i,l in enumerate(ls):
  display={'Hispanic (All Races)':'Hispanic, all races','Non-Hispanic Unknown Race':'Unknown','>=100,000':'≥100,000','2016-2019':'2016–2018'}.get(l,l).replace('2004-2009','2004–2009').replace('2010-2015','2010–2015')
  rows.append([label if i==0 else '',display]+[fmt(n(v,[l],s),s) for s in ['iCCA','eCCA']])
rad=[x['level'] for x in counts if x['variable']=='Radiation recode' and x['level'] not in ['None/Unknown','Refused (1988+)','Recommended, unknown if administered']]
for label,v,ls in [('Radiotherapy','Radiation recode',rad),('Chemotherapy','Chemotherapy recode (yes, no/unk)',['Yes']),('Registry-recorded primary sequence','Sequence number',['One primary only']),('','Sequence number',list({x['level'] for x in counts if x['variable']=='Sequence number' and x['level']!='One primary only'}))]:
 rows.append([label,'Multiple-primary sequence' if label=='' else 'Single primary' if v=='Sequence number' else 'Yes']+[fmt(n(v,ls,s),s) for s in ['iCCA','eCCA']])
assert sum(ns.values())==15578
assert sum(n('Sequence number',['One primary only'],s) for s in ['iCCA','eCCA'])==12222

(Q/'Table1_source.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
print('Saved core baseline table')
