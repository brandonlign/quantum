#!/usr/bin/env python3
"""Exploratory float screen (not proof). Pure [[n,k,d]]: iterate
  A_w forced 0 if b(n-w,d) < n+k-2w+1  (support-selection rank bound, b = LP+Gauss)
and test pure enumerator LP (MacWilliams, inclusion, shadow) for feasibility."""
import json
import numpy as np
from math import comb
from scipy.optimize import linprog
from gauss_screen import b_gauss, K
def pure_feasible(n,k,d,zero):
    M=2**(n-k); ws=[w for w in range(d,n+1) if w not in zero]
    Aeq=[[1]*len(ws)]; beq=[M-1]; Aub=[]; bub=[]
    for j in range(n+1):
        row=[K(n,j,w)/M for w in ws]; c0=K(n,j,0)/M
        sel=[1.0*(w==j) for w in ws]; s0=1.0*(j==0)
        if 1<=j<d: Aeq.append([r-s for r,s in zip(row,sel)]); beq.append(s0-c0)
        elif j>=d: Aub.append([s-r for r,s in zip(row,sel)]); bub.append(c0-s0)
        Aub.append([-K(n,j,w)*(-1)**w/M for w in ws]); bub.append(K(n,j,0)/M)
    r=linprog(np.zeros(len(ws)),A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=[(0,None)]*len(ws),method='highs')
    return r.status==0
ents=json.load(open('../../data/open_entries.json'))
res=[]
for e in ents:
    n,k,d=e['n'],e['k'],e['hi']
    if n>25: continue
    zero={w for w in range(d,n+1) if b_gauss(n-w,d) < n+k-2*w+1}
    res.append(dict(n=n,k=k,d=d,zero=sorted(zero),feasible=pure_feasible(n,k,d,zero)))
    print(res[-1],flush=True)
json.dump(res,open('forced_screen.json','w'),indent=1)
