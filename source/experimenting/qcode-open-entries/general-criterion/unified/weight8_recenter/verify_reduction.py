#!/usr/bin/env python3
"""Independent small exact checks; mathematical scope in WEIGHT8_RECENTER.md."""
import hashlib,json
from itertools import product
from math import comb
from pathlib import Path

def K(j,w):
 return sum((-1)**a*3**(j-a)*comb(w,a)*comb(9-w,j-a) for a in range(max(0,j-9+w),min(j,w)+1))
assert [[K(j,w) for w in [0,6,7,8,9]] for j in [1,2,3]] == [[27,3,-1,-5,-9],[324,-12,-12,4,36],[2268,-20,28,28,-84]]
possible=[]
for a9 in range(128):
 for a8 in range(128-a9):
  for a7 in range(128-a9-a8):
   a6=127-a7-a8-a9
   if all(K(j,0)+sum(a*K(j,w) for a,w in zip([a6,a7,a8,a9],[6,7,8,9]))>=0 for j in [1,2,3]):
    possible.append([a6,a7,a8,a9])
assert possible==[[84,0,27,16]]
assert 96 not in [0]+[sgn*2**(7-s) for s in range(4) for sgn in [-1,1]]
frames=set();count=0
for q in product(range(4),repeat=8): # 0 I,1 Z,2 X,3 Y
 a=sum(v>=2 for v in q);b=q.count(1)
 if a==0 or a%2 or a+b<6 or 8-b<6: continue
 # replacement q->q+p complements b on sites lacking X
 frames.add((a,min(b,8-a-b)));count+=1
assert frames=={(4,2),(6,0),(6,1),(8,0)}
here=Path(__file__).resolve().parent
report={'status':'PASS','integer_MacWilliams_distributions':possible,'impossible_Gauss_sum':96,'rank_two_frames':sorted(frames),'admissible_second_generators':count,'scope':'Exact finite checks supporting proof; not a target extension search or a classification audit','sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),here/'pure17_without_weight8.json',here.parent/'WEIGHT8_RECENTER.md']}}
(here/'reduction_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
