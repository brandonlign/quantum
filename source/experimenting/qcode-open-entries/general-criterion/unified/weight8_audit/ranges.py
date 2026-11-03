"""Float LP ranges of A_w for pure [[17,3,6]] under ordinary MacWilliams/inclusion/shadow (split.system, no blocks).
Exploratory only."""
import sys; sys.path.insert(0,'split')
import numpy as np
from scipy.optimize import linprog
from system import build
n,k,d=17,3,6
def rng(extra):
    w={i:0 for i in range(1,6)}; w.update(extra)
    A,b,G,cols,*_=build(n,k,d,[],weights=w)
    A=A.astype(float);b=b.astype(float);G=G.astype(float)
    out={}
    for t in range(6,18):
        c=np.array([1.0 if col[-1]==t else 0 for col in cols])
        lo=linprog(c,A_ub=G,b_ub=np.zeros(len(G)),A_eq=A,b_eq=b,bounds=[(0,None)]*len(cols),method='highs')
        hi=linprog(-c,A_ub=G,b_ub=np.zeros(len(G)),A_eq=A,b_eq=b,bounds=[(0,None)]*len(cols),method='highs')
        out[t]=(round(lo.fun,3) if lo.status==0 else lo.status, round(-hi.fun,3) if hi.status==0 else hi.status)
    return out
print('pure', rng({}))
print('mu7 (A6=0)', rng({6:0}))
print('mu8 (A6=A7=0)', rng({6:0,7:0}))
