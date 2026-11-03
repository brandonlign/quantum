"""Physical full S,D oracle verifies every added equation and inequality."""
import json,sys
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import dual,split,span,weight
from gauge_fix import supported_subspace
from local_dual_lp.system import build,letter_hist
from mixed_frames.stage3.common import atomic_json,sha

def main():
 controls=json.loads((ROOT/'rankone_recenter/verification.json').read_text())['positive_controls']+[dict(parameters=[8,2,3],S=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'])];out=[]
 for c in controls:
  n,k,d=c['parameters'];m=3 if n==5 else 4;f=n-m;S=c['S'];D=dual(S,n)
  H=[split(w,m,f)[0] for w in span(supported_subspace(D,n,range(m)))]
  C=[split(w,m,f)[1] for w in span(supported_subspace(D,n,range(m,n)))]
  N=[split(w,m,f)[1] for w in span(supported_subspace(S,n,range(m,n)))]
  actual_light=sum(weight(w,n) in (6,7) for w in span(S))
  actual_eight=sum(weight(w,n)==8 for w in span(S))
  E,b,G,V,_=build(n,k,d,m,letter_hist(H,m),Counter(weight(w,f) for w in C),Counter(weight(w,f) for w in N),light_budget=actual_light,weight8_min=actual_eight)
  counts=Counter()
  for w in span(S):
   p,s=split(w,m,f);t=(p&((1<<m)-1)).bit_count();z=((p>>m)&~p&((1<<m)-1)).bit_count();counts[((min(z,m-z-t),t),weight(s,f))]+=1
  # Variables are per orbit MEMBER counts, not sums of both members.
  x=[]
  for (z,t),w in V:
   raw=counts[((z,t),w)];x.append(raw//(1 if 2*z+t==m else 2));assert raw%(1 if 2*z+t==m else 2)==0
  assert all(v==0 for v in E@x-b)
  assert all(v<=0 for v in G@x)
  out.append(dict(parameters=c['parameters'],equations=len(b),inequalities=len(G),variables=len(V)))
 atomic_json(HERE/'verification.json',dict(status='PASS',controls=out,hashes={str(p):sha(p) for p in [HERE/'system.py',Path(__file__),ROOT.parents[1]/'pure/refined_split.py',ROOT.parents[1]/'split/system.py',ROOT/'model.py']}));print(json.dumps(out),flush=True)
if __name__=='__main__':main()
