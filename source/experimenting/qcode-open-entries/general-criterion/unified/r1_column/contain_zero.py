"""Exact Farkas: pure [[17,3,6]], fixed (A6,A7), low word u=Z^w (w=6/7) as block, plus
'no weight-8 stabilizer contains supp(u)' (contain_u=0). Valid in the all-r=1 case:
a weight-8 word containing a low word has local rank >=2, so a rank-two frame would be chosen."""
import sys,json;sys.path.insert(0,'split');sys.path.insert(0,'general-criterion/unified/frame_selection')
import numpy as np
from system import build,farkas_exact,patterns
n,k,d=17,3,6
def with_row(Meq,b,G,row,rhs):
    Meq=np.vstack([Meq,np.array([row],dtype=object)]);b=np.concatenate([b,np.array([rhs],dtype=object)])
    return Meq,b,G
out={}
for a6,a7 in [(0,1),(0,2),(1,0),(1,1),(2,0)]:
  for u in ([6] if a6 else [])+([7] if a7 else []):
    w={i:0 for i in range(1,6)};w[6]=a6;w[7]=a7
    A,b,G,cols,*_=build(n,k,d,[(u,['Z'*u])],weights=w)
    pats=patterns(u,['Z'*u])[0]
    row=[pats[y[0]][u]*int(y[1]==8-u) for y in cols]
    cert,msg=farkas_exact(*with_row(A,b,G,row,0))
    out[f'{a6}{a7}_u{u}']=dict(cert=bool(cert),msg=str(msg)[:80]);print(a6,a7,u,out[f'{a6}{a7}_u{u}'],flush=True)
json.dump(out,open('general-criterion/unified/r1_column/contain_zero.json','w'),indent=1)
