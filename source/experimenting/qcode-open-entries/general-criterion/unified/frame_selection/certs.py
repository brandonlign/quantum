#!/usr/bin/env python3
"""Exact Farkas certificates (Fraction-verified) for the frame-selection lemma, pure [[17,3,6]].
(a) A6>=3 impossible, A7>=3 impossible  (so A6,A7<=2)
(b) A8<=119 impossible                    (so A8>=120)
(c) for a weight-u stabilizer word w (u=6,7): #weight-8 stabilizer words whose support contains supp(w)
    is <= c_u  (certificate for >= c_u+1)."""
import sys, json
sys.path.insert(0,'split')
import numpy as np
from fractions import Fraction
from system import build, farkas_exact, patterns
n,k,d=17,3,6
pure={i:0 for i in range(1,6)}
def with_row(Meq,b,G,row,rhs,sense):
    """append row.x  (sense '>=' : row.x - s = rhs ; '<=' : row.x + s = rhs), slack s>=0"""
    Meq=np.vstack([np.hstack([Meq,np.zeros((Meq.shape[0],1),dtype=object)]),
                   np.array([list(row)+[-1 if sense=='>=' else 1]],dtype=object)])
    b=np.concatenate([b,np.array([rhs],dtype=object)])
    G=np.hstack([G,np.zeros((G.shape[0],1),dtype=object)])
    return Meq,b,G
out={}
def rec(tag,res):
    cert,msg=res; out[tag]=dict(cert=bool(cert),msg=str(msg)); print(tag,out[tag],flush=True)
    return bool(cert)
A,b,G,cols,*_=build(n,k,d,[],weights=pure)
for t,v in [(6,3),(7,3)]:
    row=[int(c[-1]==t) for c in cols]
    rec(f'A{t}>={v}',farkas_exact(*with_row(A,b,G,row,v,'>=')))
row=[int(c[-1]==8) for c in cols]
rec('A8<=119',farkas_exact(*with_row(A,b,G,row,119,'<=')))
for u in (6,7):
    Ab,bb,Gb,colsb,*_=build(n,k,d,[(u,['Z'*u])],weights=pure)
    pats=patterns(u,['Z'*u])[0]
    row=[pats[y[0]][u]*int(y[1]==8-u) for y in colsb]
    for c in range(0,12):
        if rec(f'contain{u}>={c}',farkas_exact(*with_row(Ab,bb,Gb,row,c,'>='))):break
json.dump(out,open('general-criterion/unified/frame_selection/certs.json','w'),indent=1)
