#!/usr/bin/env python3
"""Reproduce the exact certificate; verification uses pure/check_general.py."""
import json,sys
from math import lcm
from pathlib import Path
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here.parents[2]/'split'))
from system import build,farkas_exact
weights={w:0 for w in [1,2,3,4,5,8]}
A,b,G,columns,eq_labels,le_labels=build(17,3,6,[],weights=weights)
certificate,info=farkas_exact(A,b,G)
assert certificate,info
nu,lam=certificate
scale=lcm(*[v.denominator for v in nu+lam if v])
rows=[]
for labels,values in [(eq_labels,nu),(le_labels,lam)]:
 for label,v in zip(labels,values):
  if not v: continue
  if isinstance(label,str): row={'type':label}
  elif label[0]=='weight': row={'type':'weight','w':label[1],'value':weights[label[1]]}
  else: row={'type':label[0],'cell':list(label[1])}
  row['multiplier']=int(v*scale);rows.append(row)
res=[sum(int(v*scale)*A[i,j] for i,v in enumerate(nu))+sum(int(v*scale)*G[i,j] for i,v in enumerate(lam)) for j in range(len(columns))]
rhs=sum(int(v*scale)*b[i] for i,v in enumerate(nu))
spec={'name':'pure17_without_weight8','code':'[[17,3,6]]','description':'Pure stabilizer with A1..A5=0 and A8=0 is impossible.','blocks':[],'suffix_sites':17,'weight_assumptions':[[w,v] for w,v in weights.items()],'weight_ge_assumptions':[],'expected':{'columns':len(columns),'nu_rhs':int(rhs),'positive_residuals':sum(v>0 for v in res)},'certificate':rows}
(here/'pure17_without_weight8.json').write_text(json.dumps(spec,indent=2)+'\n')
print(info)
