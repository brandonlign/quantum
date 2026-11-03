"""Knuth random-probe estimate of the canonical Z tree size (weights + Gamma rank prune), m8,r1."""
import json,sys,random
import numpy as np
sys.path.insert(0,'..')
from model import ExtensionModel,nullspace,reduced_basis
S=json.load(open('../opus_parallel/suffix9_dimension5/suffix9_dim5.json'))['classes']
m,d=8,6
def rank_restricted(xs,G,b):
    A=nullspace(list(xs),b)
    rows=[sum((bin(a & sum(((G[i]>>0)&0) for i in [0]))) for a in A)] if False else None
    g=[sum(((bin(Gi & a).count('1'))&1)<<i for i,Gi in enumerate(G)) for a in A]
    M=[sum((bin(A[i]&g[j]).count('1')&1)<<j for j in range(len(A))) for i in range(len(A))]
    return len(reduced_basis(M))
def probe(model,rng,gamma=True):
    D=np.array(model.distances);N=len(D);idx=np.arange(N)
    xs=[];sums=[0];sizes=[0];lo=-1;est=1.0;tot=1.0;per=[1.0]
    for depth in range(m-1):
        ok=idx>lo
        for s,z in zip(sums,sizes):
            a=z+1;t=min(a,m-a);need=d-t
            if need>0: ok&=D[idx^s]>=need
        cand=idx[ok]
        if depth==m-2:
            cand=cand[(cand^sums[-1])>cand]
        if gamma and depth>=1:
            pass
        # gamma prune applies to child nodes: filter children
        if gamma and depth>=2:
            good=[c for c in cand if rank_restricted(xs+[int(c)],model.Gamma,model.b)<=2*(m-2-depth)] if len(cand)<4000 else None
            if good is None:
                sample=rng.sample(list(cand),400);frac=np.mean([rank_restricted(xs+[int(c)],model.Gamma,model.b)<=2*(m-2-depth) for c in sample])
                k=int(round(frac*len(cand)))
                good=[c for c in sample if rank_restricted(xs+[int(c)],model.Gamma,model.b)<=2*(m-2-depth)]
                if not good:return tot,per
                cnt=k
            else:cnt=len(good)
        else:
            good=list(cand);cnt=len(good)
        if cnt==0:return tot,per
        est*=cnt;tot+=est;per.append(est)
        x=int(rng.choice(good));xs.append(x)
        sums=sums+[s^x for s in sums];sizes=sizes+[z+1 for z in sizes];lo=x
    return tot,per
ci=int(sys.argv[1]) if len(sys.argv)>1 else 0
model=ExtensionModel(17,3,6,8,[255<<8],S[ci]['representative']['generators'])
rng=random.Random(5);T=[];P=[]
for k in range(int(sys.argv[2]) if len(sys.argv)>2 else 200):
    t,p=probe(model,rng);T.append(t);P.append(p+[0]*(8-len(p)))
P=np.array(P);print('class',ci,'mean tree',f'{np.mean(T):.3e}','per-depth',[f'{v:.2e}' for v in P.mean(0)])
