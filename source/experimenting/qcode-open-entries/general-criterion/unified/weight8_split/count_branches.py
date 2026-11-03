"""Exact Farkas: pure [[17,3,6]], fixed (A6,A7); minimum-word block Z^m. Also exact certs for A6>=3, A7>=3."""
import sys, json; sys.path.insert(0,'split')
import numpy as np
from system import build, farkas_exact
n,k,d=17,3,6
def run(blocks,w,ge=None):
    Meq,b,G,*_=build(n,k,d,blocks,weights=w)
    if ge is not None:  # convert last weight row to >= via slack
        Meq=np.hstack([Meq,np.zeros((Meq.shape[0],1),dtype=object)]); Meq[-1,-1]=-1
        G=np.hstack([G,np.zeros((G.shape[0],1),dtype=object)])
    return bool(farkas_exact(Meq,b,G)[0])
out={}
base={i:0 for i in range(1,6)}
out['A6>=3']=run([],{**base,6:3},ge=1); out['A7>=3']=run([],{**base,7:3},ge=1)
print(out,flush=True)
cases=[(0,0),(0,1),(0,2),(1,0),(1,1),(2,0)]
for a6,a7 in cases:
    m=6 if a6 else (7 if a7 else 8)
    w={**base,6:a6,7:a7}
    out[f'{a6}{a7}_m{m}']=run([(m,['Z'*m])],w)
    print(a6,a7,m,out[f'{a6}{a7}_m{m}'],flush=True)
json.dump(out,open('general-criterion/unified/weight8_split/count_branches.json','w'),indent=1)
