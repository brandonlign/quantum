"""Independent rational-certificate replay and physical coupling checks."""
import json
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
from screen_split_profiles import coupled_system, span, split, weight
from model import dual, reduced_basis, symp

HERE=Path(__file__).resolve().parent


def verify():
    saved=json.loads((HERE/'local_split_profile_screen.json').read_text())
    E,b,G,V=coupled_system(17,3,6,8,8)
    checked=[]
    for r in saved['records']:
        if r['status']!='EXCLUDED_EXACT':continue
        A={int(w):v for w,v in r['spectrum'].items()}
        rhs=list(b)+[16384*A.get(j,0) for j in range(10)]
        nu=[Fraction(v) for v in r['certificate']['nu']]
        lam=[Fraction(v) for v in r['certificate']['lam']]
        assert len(nu)==len(E) and len(lam)==len(G) and all(x>=0 for x in lam)
        en=[(i,x) for i,x in enumerate(nu) if x]
        gn=[(i,x) for i,x in enumerate(lam) if x]
        coefficients=[sum(x*int(E[i,j]) for i,x in en)+sum(x*int(G[i,j]) for i,x in gn)
                      for j in range(len(V))]
        value=sum(nu[i]*int(rhs[i]) for i in range(len(rhs)))
        assert all(c>=0 for c in coefficients) and value<0
        checked.append(dict(class_indices=r['class_indices'],rhs=str(value),
                            minimum_coefficient=str(min(coefficients))))

    # Check the added MacWilliams equalities on actual target-size isotropic
    # spaces. They need not meet the target distance: the coupling alone is
    # a structural identity, independent of purity/minimum assumptions.
    rng=random.Random(1736);physical=[];n=17;m=8;f=9;p=((1<<m)-1)<<n
    for trial in range(8):
        S=[1<<(n+i) for i in range(14)]
        for _ in range(70):
            v=rng.randrange(1<<(2*n))
            if symp(v,p,n):v^=1
            S=[s^(v if symp(s,v,n) else 0) for s in S]
        assert len(reduced_basis(S))==14
        assert all(symp(a,b,n)==0 for a in S for b in S)
        assert p in span(S)
        P=reduced_basis([split(s,m,f)[1] for s in S])
        C=dual(P,f);Aw=Counter(weight(c,f) for c in span(C))
        counts=Counter()
        for s in span(S):
            prefix,suffix=split(s,m,f);x=prefix&((1<<m)-1);z=prefix>>m
            counts[((z&~x).bit_count(),x.bit_count(),weight(suffix,f))]+=1
        variables=np.array([counts[(z,t,w)] for ((z,t),w) in V],dtype=object)
        expected=np.array([16384*Aw.get(j,0) for j in range(10)],dtype=object)
        assert all(v==0 for v in E[-10:]@variables-expected)
        physical.append(dict(trial=trial,local_rank=14-len(P),suffix_dimension=len(C),
                             stabilizer_words=16384))
    result=dict(status='PASS',certificate_replays=checked,physical_coupling_checks=physical,
                scope='Exact rational replay without system.verify and physical MacWilliams coupling at target dimensions')
    (HERE/'split_profile_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':verify()
