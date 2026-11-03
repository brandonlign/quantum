"""Physical verification of minimum-support projection-rank bounds."""
import hashlib
import json
import random
from itertools import combinations
from pathlib import Path
from collections import Counter
from model import span,weight,reduced_basis,dual,symp,embed
from projection_rank_bound import projection_rank_bound


def restrict(w,sites,n):
    return sum(((w>>i)&1)<<j for j,i in enumerate(sites)) | (
        sum(((w>>(n+i))&1)<<j for j,i in enumerate(sites))<<len(sites))


def check(S,n,T=None):
    S=reduced_basis(S);words=span(S)
    assert all(symp(a,b,n)==0 for a in S for b in S)
    m=min(weight(w,n) for w in words if w)
    if T is None:
        p=next(w for w in words if w and weight(w,n)==m)
        support=(p&((1<<n)-1))|(p>>n)
        T=[i for i in range(n) if support>>i&1]
    assert len(T)==m
    W=[i for i in range(n) if i not in T];f=len(W)
    R=reduced_basis([w for w in words if restrict(w,W,n)==0]);r=len(R)
    assert r in [1,2] and all(weight(w,n)==m for w in span(R) if w)
    P=reduced_basis([restrict(s,W,n) for s in S]);C=dual(P,f)
    checks=nontrivial=0
    for t in range(min(f,3)+1):
        for A in combinations(range(f),t):
            actual=len(reduced_basis([restrict(c,A,f) for c in C]))
            bound=projection_rank_bound(m,r,t)
            assert actual>=bound['required_projection_rank'],(n,S,T,A,actual,bound)
            checks+=1;nontrivial+=bound['required_projection_rank']>0
    return dict(n=n,m=m,r=r,suffix_dimension=len(C),projection_checks=checks,
                nontrivial_checks=nontrivial,suffix_generators=C)


def verify():
    rng=random.Random(131736);counts=Counter();checks=nontrivial=0
    for _ in range(1000):
        n=rng.randrange(4,11);s=rng.randrange(1,min(n,5)+1)
        S=[1<<(n+i) for i in range(s)]
        for _ in range(35):
            v=rng.randrange(1<<(2*n))
            S=[w^(v if symp(w,v,n) else 0) for w in S]
        r=check(S,n);counts[(r['m'],r['r'])]+=1
        checks+=r['projection_checks'];nontrivial+=r['nontrivial_checks']
    # Sharp rank-one even boundary: minimum-six stabilizer on seven sites,
    # with a rank-one six-site local stabilizer and zero suffix code.
    even6=check([8064,2127,8892],7,list(range(6)))
    assert even6['m']==6 and even6['r']==1 and even6['suffix_dimension']==0
    assert projection_rank_bound(6,1,1)['required_projection_rank']==0
    # Odd minimum-five rank-one case: rank-one suffix projection attains
    # the necessary lower bound, while zero projection is forbidden.
    odd5=check([1984,47],6,list(range(5)))
    assert odd5['m']==5 and odd5['r']==1 and odd5['suffix_dimension']==1
    assert projection_rank_bound(5,1,1)['required_projection_rank']==1
    rank2=check([63,16128,15567],8,list(range(6)))
    assert rank2['m']==6 and rank2['r']==2
    assert projection_rank_bound(6,2,1)['required_projection_rank']==2
    def pauli(text):
        n=len(text)
        return sum((c in 'XY')<<i for i,c in enumerate(text)) | (sum((c in 'ZY')<<i for i,c in enumerate(text))<<n)
    perfect5=check([pauli(t) for t in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']],5,list(range(4)))
    assert perfect5['m']==4 and perfect5['r']==2 and perfect5['suffix_dimension']==0
    assert projection_rank_bound(4,2,1)['required_projection_rank']==0
    source=Path(__file__).resolve().parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    classes=json.loads(source.read_text())['classes'];excluded=[]
    for c in classes:
        C=c['representative']['generators'];support=0
        for w in C:support|=(w&511)|(w>>9)
        if support!=511:
            assert projection_rank_bound(8,1,1)['required_projection_rank']==1
            excluded.append(dict(cls=c['class_index'],zero_sites=[i for i in range(9) if not (support>>i&1)]))
    assert excluded==[dict(cls=24383,zero_sites=[7])]
    result=dict(status='PASS',physical_codes=1000,projection_checks=checks,nontrivial_checks=nontrivial,
                minimum_support_types=[dict(m=m,r=r,count=count) for (m,r),count in sorted(counts.items())],
                controls=dict(even_rankone_boundary=even6,odd_rankone_sharp=odd5,
                              ranktwo_full_projection=rank2,perfect_five_boundary=perfect5),
                target_analytic_exclusions=excluded,
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [source,Path(__file__),Path(__file__).with_name('projection_rank_bound.py')]},
                scope='Necessary classical projection bounds on physical isotropic stabilizers; target exclusion requires no extension search')
    Path(__file__).with_name('projection_rank_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k in ['status','physical_codes','projection_checks','nontrivial_checks','target_analytic_exclusions']})


if __name__=='__main__':verify()
