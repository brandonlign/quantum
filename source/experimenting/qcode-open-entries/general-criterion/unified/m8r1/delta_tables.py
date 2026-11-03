"""Coset-distance tables delta(s) for every (9,2^5,>=6) suffix class (X|Z layout).

Syndrome s(w)_j=<h_j,w> for a fixed basis h_0..h_12 of P=C^perp (model.dual).
Writes delta_tables.npy (uint8, classes x 8192) and P bases to delta_meta.json.
"""
import json,sys,time,hashlib
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from model import dual,symp,reduced_basis,span
f=9
src=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
S=json.loads(src.read_text())
errs=[]
for i in range(f):
    for x,z in [(1,0),(0,1),(1,1)]:errs.append((x<<i)|(z<<(f+i)))
T=np.zeros((len(S['classes']),1<<13),np.uint8);bases=[];t0=time.time()
for ci,c in enumerate(S['classes']):
    C=c['representative']['generators'];assert len(reduced_basis(C))==5
    P=dual(C,f);assert len(P)==13 and all(symp(p,g,f)==0 for p in P for g in C)
    syn=np.array([sum(symp(h,e,f)<<j for j,h in enumerate(P)) for e in errs],np.int64)
    d=np.full(1<<13,255,np.uint8);d[0]=0;fr=np.array([0]);lvl=0
    while fr.size:
        lvl+=1;nb=np.unique((fr[:,None]^syn[None,:]).ravel());nb=nb[d[nb]==255];d[nb]=lvl;fr=nb
    assert (d<255).all();T[ci]=d;bases.append(P)
    if ci%5000==0:print(ci,round(time.time()-t0,1),flush=True)
np.save(HERE/'delta_tables.npy',T)
(HERE/'delta_meta.json').write_text(json.dumps(dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),layout='x|z<<9',P_bases=bases))+'\n')
deep=(T>=5).sum(1);print('deep>=5 count distribution',np.bincount(deep)[:40].tolist(),'max',deep.max())
print('covering radius distribution',np.bincount(T.max(1)).tolist())
