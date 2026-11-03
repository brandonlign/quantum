"""Root-session independent check of the mixed-frame capacity prefilter (numpy brute force)."""
import json, numpy as np
d=json.load(open('suffix9_classes.json'))['classes']
W=np.arange(1<<18,dtype=np.int64)
pc=np.vectorize(lambda v:bin(v).count('1'))
wt=np.array([bin((w&511)|(w>>9)).count('1') for w in range(1<<18)],dtype=np.int64)
def sp_all(g):  # symplectic product of every word with g
    a=(W&511)&(g>>9); b=(W>>9)&(g&511)
    return (pc8(a)+pc8(b))&1
POP=np.array([bin(i).count('1') for i in range(512)],dtype=np.int64)
def pc8(v): return POP[v]
out={}
for c in d:
    gens=c['representative']['generators']
    mask=np.ones(1<<18,bool)
    for g in gens: mask&=(sp_all(g)==0)
    perp=W[mask]; assert len(perp)==1<<12
    basis=[]
    for v in perp.tolist():
        r=v
        for b in basis: r=min(r,r^b)
        if r: basis.append(r)
        if len(basis)==12: break
    syn=np.zeros(1<<18,dtype=np.int64)
    for j,b in enumerate(basis): syn|=sp_all(b)<<j
    delta=np.full(4096,99); np.minimum.at(delta,syn,wt)
    out[int(c['class_index'])]=int((delta>=5).sum())
rej={'668':sorted(k for k,v in out.items() if v<4),'688':sorted(k for k,v in out.items() if v<2)}
print(rej); json.dump(dict(counts=out,rejects=rej),open('mixed_frames_root_check.json','w'),indent=1)
