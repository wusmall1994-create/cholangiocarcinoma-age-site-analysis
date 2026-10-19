import os,gzip
from pathlib import Path
import csv,json,math
from collections import defaultdict
Q=Path(os.environ['CCA_RESULTS_DIR'])
groups=defaultdict(list); byrep=defaultdict(dict)
matrix_path=Q/'aggregate_bootstrap_contrasts.csv'
if not matrix_path.exists(): matrix_path=Q/'aggregate_bootstrap_contrasts.csv.gz'
open_matrix=lambda: gzip.open(matrix_path,'rt',encoding='utf-8-sig',newline='') if matrix_path.suffix=='.gz' else matrix_path.open(encoding='utf-8-sig',newline='')
with open_matrix() as f:
 for r in csv.DictReader(f):
  k=tuple(r[x] for x in ['sequence_group','definition','state','contrast']);v=float(r['estimate']);b=int(r['replicate'])
  assert math.isfinite(v);groups[k].append(v);byrep[b][k]=v
def quantile7(a,p):
 a=sorted(a);h=(len(a)-1)*p;i=math.floor(h);return a[i]+(h-i)*(a[min(i+1,len(a)-1)]-a[i])
errors=[];identity=[]
with (Q/'ridge_contrasts_with_CI.csv').open(encoding='utf-8-sig',newline='') as f:
 for r in csv.DictReader(f):
  k=tuple(r[x] for x in ['sequence_group','definition','state','contrast']);a=groups[k];assert len(a)==2000
  for field,p in [('lower95',.025),('upper95',.975)]:errors.append(abs(float(r[field])-quantile7(a,p)))
for b,r in byrep.items():
 assert len(r)==len(groups)
 for (g,d,s,c),v in r.items():
  if c=='RD_change_75_minus_50':
   identity.append(abs(v-(r[(g,d,s,'RD_age75')]-r[(g,d,s,'RD_age50')])))
  if d.endswith('_minus_SEER'):
   base=d.removesuffix('_minus_SEER');identity.append(abs(v-(r[(g,base,s,c)]-r[(g,'SEER_original',s,c)])))
  if s=='all_death' and d!='Direct_allcause_including_unknown':
   identity.append(abs(v-(r[(g,d,'cancer',c)]+r[(g,d,'other',c)])))
assert len(byrep)==2000 and max(errors)<1e-12 and max(identity)<1e-12
report={'status':'PASS','replicates':len(byrep),'contrasts_per_replicate':len(groups),'CI_endpoints_checked':len(errors),'max_CI_error':max(errors),'identities_checked':len(identity),'max_identity_error':max(identity),'method':'Independent Python CSV parsing, manual type-7 interpolation; paired endpoint, age and complementary-state identities'}
(Q/'independent_bootstrap_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report))
