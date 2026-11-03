#!/usr/bin/env python3
"""Split relaxation strengthened by exact local data for a weight-eight frame of pure [[17,3,6]]:
  stabilizers supported on T:  sum_{c_f=0} A_c = 2^r     (S cap V_T = R exactly)
  centralizer supported on W:  sum_{c_T=0} B_c = 2^(4+r) (dim C = n+k-2m+r)
CERT = exactly verified infeasibility; 'none' proves nothing."""
import sys, json, itertools
sys.path.insert(0, 'split')
import numpy as np
from system import build, farkas_exact, kraw
n,k,d=17,3,6; M=2**(n-k)
frames={'r1':['ZZZZZZZZ'],'r2_668':['ZZZZZZZZ','XXXXZZII'],'r2_688':['ZZZZZZZZ','XXXXXXII'],
        'r2_778':['ZZZZZZZZ','XXXXXXZI'],'r2_888':['ZZZZZZZZ','XXXXXXXX']}
def strengthened(g, extra_w=None):
    r=len(g); w={i:0 for i in range(1,6)}
    if extra_w: w.update(extra_w)
    Meq,b,G,cols,le,ll=build(n,k,d,[(8,g)],weights=w)
    from system import patterns
    pats=patterns(8,g)[0]; f=9
    Tb=[[sum(kraw(8,a,j)*P[j] for j in range(9)) for a in range(9)] for P in pats]
    rowA=[ (pats[y[0]][0]*0 + sum(pats[y[0]][a] for a in range(9))) * (y[1]==0) for y in cols]
    rowB=[ Tb[y[0]][0]*kraw(f,0,y[1]) * 1 for y in cols]  # sum over suffix cells c_f of T(0,c_f;y): handled below
    rowB=[sum(Tb[y[0]][0]*kraw(f,cf,y[1]) for cf in range(f+1)) for y in cols]
    Meq=np.vstack([Meq,np.array([rowA],dtype=object),np.array([rowB],dtype=object)])
    b=np.concatenate([b,np.array([2**r, M*2**(4+r)],dtype=object)])
    return Meq,b,G
out={}
for tag,g in frames.items():
    Meq,b,G=strengthened(g)
    cert,msg=farkas_exact(Meq,b,G)
    out[tag]=dict(cert=bool(cert),msg=str(msg)); print(tag,out[tag],flush=True)
json.dump(out,open('general-criterion/unified/weight8_split/local_exact.json','w'),indent=1)
