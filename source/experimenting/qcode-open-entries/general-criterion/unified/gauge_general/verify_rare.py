"""Physical singleton basis and complete S6 x tail-S2 normal-form coverage."""
import itertools,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import ExtensionModel,coordinates,lift_coordinate,span,symp,weight
from gauge_general.rare_frame import rare_frame
from gauge_general.common import quotient
from gauge_general.verify_backend import native
from mixed_frames.stage3.common import atomic_json,sha

def permute(w,sites):return sum(((w>>j)&1)<<sites[j] for j in range(8))|sum(((w>>(8+j))&1)<<(8+sites[j]) for j in range(8))
def canonical(xs,body):
 tails=8-body
 if body==6 and xs[1]<=xs[0]:return False
 for depth in range(tails+1,7):
  cur=list(xs[tails:depth])
  if cur!=sorted(cur):return False
  if depth==6:
   last=0
   for x in xs:last^=x
   if last<=cur[-1]:return False
   cur.append(last)
  points=[0]+cur
  if any(sorted(v^a for v in points if v!=a)<cur for a in cur):return False
 return True

def main():
 report=[]
 for name in ('668','688','778'):
  f,body=rare_frame(name);physical=f.verify_physical();fixture=[1<<i for i in range(6)];actions=set()
  # Body permutations preserve every physical R. Tail swap for 778 is
  # accompanied by phase on each body site, preserving p and swapping q/p+q.
  tails=[False,True] if body==6 else [False]
  for order in itertools.permutations(range(body)):
   for swap in tails:
    sites=list(order)+list(range(body,8))
    if swap:sites[-2:]=sites[-2:][::-1]
    def apply(w):
     v=permute(w,sites)
     if name=='778' and swap:
      v^=(v&((1<<body)-1))<<8
     return v
    assert {apply(w) for w in span(f.R)}==set(span(f.R))
    units=[1<<i for i in range(16)];assert all(symp(apply(a),apply(b),8)==symp(a,b,8) for a in units for b in units)
    A=[coordinates(apply(w),f.R+f.E)>>2 for w in f.E[:6]];assert all(a<64 for a in A)
    actions.add(tuple(lift_coordinate(a,fixture) for a in A))
  expected=24 if body==4 else 1440;assert len(actions)==expected
  accepted=[xs for xs in actions if canonical(xs,body)]
  assert len(accepted)==1
  model=ExtensionModel(14,0,1,8,f.R,[]);data=quotient(model)
  samples=accepted+sorted(set(actions)-set(accepted))[:12]
  for xs in samples:
   st,leaves=native(model,f,data,fast=True,**{'prefix':','.join(map(str,xs)),'body-circuit':body})
   assert st['status']=='COMPLETE' and bool(leaves)==canonical(xs,body)
  report.append(dict(frame=name,physical_words=physical,actions=len(actions),representatives=len(accepted),native_samples=len(samples)))
 atomic_json(HERE/'rare_verification.json',dict(status='PASS',frames=report,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'rare_frame.py',HERE/'enumerate.cpp',HERE/'enumerate',ROOT/'mixed_frames/frames.py',ROOT/'model.py']}));print(json.dumps(report),flush=True)
if __name__=='__main__':main()
