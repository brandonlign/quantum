#!/usr/bin/env python3
"""Exploratory float screen. For pure [[n,k,d]]: a stabilizer word of weight w forces an
additive (n-w, 2^(n+k-2w+r), >=d) code with r>=1. If b(n-w,d) < n+k-2w+1 then A_w=0.
Add those zeros to the pure enumerator LP. Nothing here is a proof (float LP)."""
import json, math, sys
from math import comb
import numpy as np
from scipy.optimize import linprog
def K(n,j,w,q=4):
    return sum((-1)**a*(q-1)**(j-a)*comb(w,a)*comb(n-w,j-a) for a in range(max(0,j-n+w),min(j,w)+1))
def add_lp(f,d):
    """LP upper bound on |C| for additive GF(4) length f distance d (A_0=1)."""
    if f<d: return 1.0
    ws=list(range(d,f+1))
    c=-np.ones(len(ws))
    A=[[-K(f,j,w) for w in ws] for j in range(1,f+1)]; b=[K(f,j,0) for j in range(1,f+1)]
    r=linprog(c,A_ub=A,b_ub=b,bounds=[(0,None)]*len(ws),method='highs')
    return 1-r.fun
def bdim(f,d): return math.floor(math.log2(add_lp(f,d))+1e-9)
def pure_lp(n,k,d,zero):
    """feasibility of pure stabilizer enumerator: A_0=1, A_w=0 for w<d and w in zero; B=MacW>=0, B_w=A_w for w<d,
    B_w>=A_w, shadow>=0. Variables A_d..A_n."""
    M=2**(n-k); ws=[w for w in range(d,n+1) if w not in zero]
    if not ws: return False
    Aeq=[[1]*len(ws)]; beq=[M-1]; Aub=[]; bub=[]
    for j in range(0,n+1):
        row=[K(n,j,w)/M for w in ws]; const=K(n,j,0)/M   # B_j = (1/M)sum A_w K_j(w) * 2^? 
        # normalizer enumerator B_j = 2^{-(n-k)} * ... ; |S^perp|=2^{n+k}: B_j=(2^{n+k}/4^n)... use B_j = (1/2^{n-k}) sum_w A_w K_j(w)
        Aj=[1 if (w==j) else 0 for w in ws]; aj=1 if j==0 else 0
        if j<d:
            if j>0: Aeq.append([r-a for r,a in zip(row,Aj)]); beq.append(aj-const)
        else:
            Aub.append([a-r for r,a in zip(row,Aj)]); bub.append(const-aj)
        srow=[K(n,j,w)*(-1)**w/M for w in ws]; Aub.append([-s for s in srow]); bub.append(K(n,j,0)/M)
    # scale
    r=linprog(np.zeros(len(ws)),A_ub=np.array(Aub),b_ub=bub,A_eq=np.array(Aeq),b_eq=beq,bounds=[(0,None)]*len(ws),method='highs')
    return r.status==0
ents=json.load(open('../../data/open_entries.json'))
out=[]
for e in ents:
    n,k,d=e['n'],e['k'],e['hi']
    if n>30: continue
    base=pure_lp(n,k,d,set())
    zero={w for w in range(d,n+1) if bdim(n-w,d) < n+k-2*w+1}
    new=pure_lp(n,k,d,zero) if base else None
    out.append(dict(n=n,k=k,d=d,pure_lp_feasible=base,forced_zero=sorted(zero),with_zero=new))
    if base and not new: print('CANDIDATE',n,k,d,sorted(zero),flush=True)
json.dump(out,open('explore.json','w'),indent=1)
print(sum(o['pure_lp_feasible'] for o in out),'pure-feasible of',len(out))
