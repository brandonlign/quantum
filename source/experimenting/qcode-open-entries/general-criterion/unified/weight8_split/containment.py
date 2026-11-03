"""Float LP (exploratory): for pure [[17,3,6]] with fixed (A6,A7), a low word w (weight u) on block U.
max # weight-8 stabilizers whose support contains U  = A_{(u, 8-u)} in split cells; and min A8."""
import sys; sys.path.insert(0,'split')
import numpy as np
from scipy.optimize import linprog
from system import build
n,k,d=17,3,6
def lp(blocks,w,target):
    Meq,b,G,cols,*_=build(n,k,d,blocks,weights=w)
    from system import patterns
    Mf,bf,Gf=Meq.astype(float),b.astype(float),G.astype(float)
    c=np.array([target(y) for y in cols],float)
    lo=linprog(c,A_ub=Gf,b_ub=np.zeros(len(Gf)),A_eq=Mf,b_eq=bf,bounds=[(0,None)]*len(cols),method='highs')
    hi=linprog(-c,A_ub=Gf,b_ub=np.zeros(len(Gf)),A_eq=Mf,b_eq=bf,bounds=[(0,None)]*len(cols),method='highs')
    return (lo.fun if lo.status==0 else None, -hi.fun if hi.status==0 else None)
for a6,a7 in [(0,1),(0,2),(1,0),(1,1),(2,0)]:
    w={i:0 for i in range(1,6)}; w[6]=a6; w[7]=a7
    A8=lp([],w,lambda y: 1.0*(y[-1]==8))
    res={}
    for u in ([6] if a6 else [])+([7] if a7 else []):
        from system import patterns
        pats=patterns(u,['Z'*u])[0]
        res[u]=lp([(u,['Z'*u])],w,lambda y,u=u,pats=pats: pats[y[0]][u]*(y[1]==8-u))
    print((a6,a7),'A8 range',A8,'max containing',res,flush=True)
