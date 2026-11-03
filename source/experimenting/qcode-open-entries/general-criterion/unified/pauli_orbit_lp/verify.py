"""Physical controls check the full local orbit transform independently."""
import json,sys
from pathlib import Path
from collections import Counter
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import dual,span,split,weight,coordinates
from gauge_fix import supported_subspace
from pauli_orbit_lp.system import prepare,build
from mixed_frames.stage3.common import atomic_json,sha

def main():
 saved=json.loads((ROOT/'rankone_recenter/verification.json').read_text())['positive_controls'][1]
 c8=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'];records=[]
 for n,k,S in [(7,1,saved['S']),(8,2,c8)]:
  m=4;f=n-m;D=dual(S,n)
  H=[split(w,m,f)[0] for w in supported_subspace(D,n,range(m))];C=[split(w,m,f)[1] for w in supported_subspace(D,n,range(m,n))];N=[split(w,m,f)[1] for w in supported_subspace(S,n,range(m,n))]
  data=prepare(H,m);E,b,G,_=build(data,n,k,3,Counter(weight(w,f) for w in span(C)),Counter(weight(w,f) for w in span(N)))
  counts=Counter()
  for w in span(S):
   p,s=split(w,m,f);a=coordinates(p,[data['p']]+data['A'])>>1;counts[(a,weight(s,f))]+=1
  # Average counts for one physical prefix member over its A-orbit.
  x=np.array([sum(counts[(a,b)] for a in o)/(2*len(o)) for o in data['AO'] for b in range(f+1)])
  assert np.max(np.abs(E@x-b))<1e-9 and np.max(G@x)<=1e-9
  # Independent direct D suffix counts at every B representative, then
  # average over B-orbits to compare the exact character transform.
  dc=Counter()
  for w in span(D):
   p,s=split(w,m,f);index=sum(__import__('model').symp(a,p,m)<<i for i,a in enumerate(data['A']));dc[(index,weight(s,f))]+=1
  from system import kraw
  K=np.array([[kraw(f,j,b) for b in range(f+1)] for j in range(f+1)],dtype=np.int64);T=np.kron(data['char'],K);actual=T@x/(1<<(n-k-1))
  expected=np.array([sum(dc[(v,j)] for v in o)/(len(o)*(1<<len(H))) for o in data['BO'] for j in range(f+1)])
  assert np.max(np.abs(actual-expected))<1e-9
  records.append(dict(parameters=[n,k,3],A_orbits=len(data['AO']),B_orbits=len(data['BO']),physical_order=data['physical_order'],equations=len(b),inequalities=len(G)))
 atomic_json(HERE/'verification.json',dict(status='PASS',controls=records,hashes={str(p):sha(p) for p in [Path(__file__),HERE/'system.py',ROOT/'model.py',ROOT/'m8r1/suffix_symmetries.py',ROOT.parents[1]/'split/system.py']}));print(json.dumps(records),flush=True)
if __name__=='__main__':main()
