"""Complete independent label oracle plus pure non-state positive controls."""
import itertools,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import embed,reduced_basis,span,weight,lift_coordinate,dual
from gauge_fix import fix_prefix_kernel
from subcode_selection import SelectionContext
from gauge_general.inverse import ordered_basis,select
from mixed_frames.stage3.common import atomic_json,sha

def brute(S,n,m,d,cut,budget):
 B,h=ordered_basis(S,n,m);ctx=SelectionContext(B,n)
 assert min(weight(w,n) for w in span(dual(B,n)) if w)>=d
 domains=[[v for v in range(1<<len(B)) if v&((1<<h)-1)==1<<(j+1)] for j in range(h-1)]
 low=[c for c,w in enumerate(ctx.weights) if c and w<cut];out=[]
 for rows in itertools.product(*domains):
  if any(ctx.distances[v]<d for v in span(rows) if v):continue
  if sum(not any((c&v).bit_count()%2 for v in rows) for c in low)>budget:continue
  out.append(list(rows))
 return out

def main():
 saved=json.loads((ROOT/'rankone_recenter/verification.json').read_text());S7=saved['positive_controls'][1]['S'];fixed=fix_prefix_kernel(S7,7,range(4));fixtures=[(fixed['S_fixed'],7,4,3)]
 for m,h,n in [(4,2,6),(4,3,6),(5,4,7)]:
  f=n-m;R=[((1<<m)-1)<<m]+[((1<<i)|(1<<(i+1)))<<m for i in range(h-1)]
  S=[embed(v,0,m,f) for v in R]+[embed(0,1<<f,m,f)]
  fixtures.append((S,n,m,1))
 records=[]
 for S,n,m,d in fixtures:
  for cut,budget in [(d+1,0),(d+1,2),(d+1,1000),(n+1,0),(n+1,1000)]:
   expected=brute(S,n,m,d,cut,budget);native=select(S,n,m,d,cut,budget)
   assert (native['status']=='FOUND')==bool(expected)
   if expected:assert native['rows'] in expected
   else:assert native['status']=='EXCLUDED_BASE'
   records.append(dict(n=n,rank=len(S),h=ordered_basis(S,n,m)[1],d=d,cut=cut,budget=budget,admissible=len(expected),status=native['status']))
  assert select(S,n,m,d,n+1,0,nodes=0)['status']=='LIMIT'
 data=json.loads((ROOT/'radical_projection/positive_8.json').read_text());original=data['result']['generators'];state=data['state_basis']
 base=reduced_basis(original+[state[1]]);assert len(base)==7
 answer=select(base,8,4,3,4,1000);assert answer['status']=='FOUND'
 # Tensor with the standard pure five-qubit [[5,1,3]] stabilizer.
 letters={'I':0,'X':1,'Z':2,'Y':3}
 def word(s):return sum((letters[v]&1)<<j for j,v in enumerate(s))|(sum((letters[v]>>1)<<j for j,v in enumerate(s))<<len(s))
 S5=[word(s) for s in ('XZZXI','IXZZX','XIXZZ','ZXIXZ')]
 assert min(weight(w,5) for w in span(dual(S5,5)) if w)>=3
 base13=[embed(w,0,8,5) for w in state]+[embed(0,w,8,5) for w in S5]
 a13=select(base13,13,4,3,4,10000);assert a13['status']=='FOUND' and len(a13['generators'])==10
 report=dict(status='PASS',complete_fixtures=records,pure_nonstate_controls=[[8,2,3],[13,3,3]],hashes={str(p.relative_to(ROOT)):sha(p) for p in [HERE/'inverse.cpp',HERE/'inverse',HERE/'inverse.py',Path(__file__),ROOT/'model.py',ROOT/'subcode_selection.py']})
 atomic_json(HERE/'verification.json',report);print(json.dumps(dict(status='PASS',complete_fixtures=len(records),pure_nonstate_controls=report['pure_nonstate_controls'])),flush=True)
if __name__=='__main__':main()
