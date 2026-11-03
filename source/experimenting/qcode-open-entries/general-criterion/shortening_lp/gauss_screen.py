#!/usr/bin/env python3
"""Exploratory float screen (not proof).
b(f,d): largest c such that an additive (f,2^c,>=d) distribution passes MacWilliams LP together with
the Gauss-sum condition sum(-1)^w A_w in {0,+-2^(c-s)} (some s). Then for each open [[n,k,d]], pure
case: find weights w forced (LP over pure enumerators) to have A_w>0 or allowed, and local-rank bound
r <= b(n-w,d)-n-k+2w. Report w with bound<1 (then A_w=0 forced) and re-run LP."""
import json, math
from math import comb
from functools import lru_cache
import numpy as np
from scipy.optimize import linprog
def K(n,j,w): return sum((-1)**a*3**(j-a)*comb(w,a)*comb(n-w,j-a) for a in range(max(0,j-n+w),min(j,w)+1))
def feas(f,d,c,S=None):
    ws=list(range(d,f+1))
    if not ws: return c==0
    Aeq=[[1]*len(ws)]; beq=[2**c-1]
    if S is not None: Aeq.append([(-1)**w for w in ws]); beq.append(S-1)
    Aub=[[-K(f,j,w) for w in ws] for j in range(1,f+1)]; bub=[K(f,j,0) for j in range(1,f+1)]
    r=linprog(np.zeros(len(ws)),A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=[(0,None)]*len(ws),method='highs')
    return r.status==0
@lru_cache(None)
def b_lp(f,d):
    c=0
    while c<=2*f and feas(f,d,c+1): c+=1
    return c
@lru_cache(None)
def b_gauss(f,d):
    c=b_lp(f,d)
    while c>0:
        Ss={0}|{sg*2**(c-s) for s in range(c//2+1) for sg in (1,-1)}
        if any(feas(f,d,c,S) for S in Ss): return c
        c-=1
    return 0
diff=[(f,d,b_lp(f,d),b_gauss(f,d)) for d in range(3,9) for f in range(d,22) if b_lp(f,d)!=b_gauss(f,d)]
print('Gauss improvements (f,d,LP,LP+Gauss):',diff)
json.dump(dict(improvements=diff),open('gauss_bounds.json','w'),indent=1)
