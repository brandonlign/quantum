"""Direct restricted-Gram check of the partial rank-budget identity."""
import hashlib
import json
import random
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from model import reduced_basis,nullspace,lift_coordinate


def verify():
    rng=random.Random(1736);checks=0
    for b in range(2,14):
        for _ in range(150):
            G=[0]*b
            for i in range(b):
                for j in range(i):
                    if rng.randrange(2):G[i]^=1<<j;G[j]^=1<<i
            def pair(u,v):
                return (u&sum(((r&v).bit_count()%2)<<i for i,r in enumerate(G))).bit_count()%2
            av=[1<<i for i in range(b)];U=[]
            while av:
                ij=next(((i,j) for i in range(len(av)) for j in range(i+1,len(av))
                         if pair(av[i],av[j])),None)
                if ij is None:break
                i,j=ij;u,v=av[i],av[j];U.extend([u,v])
                av=[w^(u if pair(w,v) else 0)^(v if pair(w,u) else 0)
                    for k,w in enumerate(av) if k not in (i,j)]
            R=av;rho=len(U);h=len(R)
            X=reduced_basis([rng.randrange(1<<b) for _ in range(rng.randrange(b+1))]);z=len(X)
            rx=[sum(((x&r).bit_count()%2)<<i for i,r in enumerate(R)) for x in X]
            u=len(reduced_basis(rx))
            # Independent coefficient-space nullspace, rather than the C
            # implementation's incremental elimination with carried Q-values.
            rel=nullspace([sum(((v>>i)&1)<<j for j,v in enumerate(rx)) for i in range(h)],z)
            zero=[lift_coordinate(c,X) for c in rel]
            Q=[sum(((x&v).bit_count()%2)<<i for i,v in enumerate(U)) for x in zero]
            def qp(v,w):
                return sum(((v>>i)&1)*((w>>(i^1))&1) for i in range(rho))%2
            F=[sum(qp(v,w)<<j for j,w in enumerate(Q)) for v in Q]
            expected=rho-2*z+2*u+len(reduced_basis(F))
            A=nullspace(X,b)
            actual=len(reduced_basis([sum(pair(v,w)<<j for j,w in enumerate(A)) for v in A]))
            assert expected==actual,(b,X,G,expected,actual)
            checks+=1
    p=Path(__file__).resolve()
    result=dict(status='PASS',partial_rank_identities=checks,
                sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                scope='Direct restricted Gram rank versus radical-evaluation kernel formula')
    p.with_name('rank_identity_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':verify()
