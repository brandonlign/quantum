"""Finite checks supporting the radical-pair spectrum lemma."""
from collections import Counter
import hashlib
import json
from math import comb
from pathlib import Path
from model import symp, weight, span, reduced_basis, gram

BASE = Path(__file__).resolve().parent


def verify():
    compositions = 0
    for s in [2, 4, 6, 8]:
        m = 4 * s
        nlocal=m+1
        kraw=lambda j,w:sum((-1)**i*3**(j-i)*comb(w,i)*comb(nlocal-w,j-i)
                           for i in range(max(0,j-(nlocal-w)),min(j,w)+1))
        for w in range(nlocal+1):
            assert 8*(w-3*s)*(w-3*s-2)==kraw(2,w)+2*kraw(1,w)+6*(s-1)
        allowed_profiles={(3*s,)*4,(3*s+1,)*4,(3*s,3*s,3*s+2,3*s+2)}
        if s>=4:
            allowed_profiles.add((3*s,3*s,3*s,3*s+4))
        for a in range(m+1):
            for b in range(m-a+1):
                for c in range(m-a-b+1):
                    counts = [a,b,c,m-a-b-c]
                    if len({v%2 for v in counts}) != 1:
                        continue
                    for tail in [0,1]:
                        weights = [m-v+tail for v in counts]
                        if min(weights)<3*s:
                            continue
                        compositions += 1
                        assert tuple(sorted(weights)) in allowed_profiles
                        if weights[0]%2:
                            assert weights == [3*s+1]*4
    n=9
    R=[255,255<<n]
    shifts=[0,R[0],R[1],R[0]^R[1]]
    deep, odd = Counter(), Counter()
    cosets=0
    for w in range(1<<(2*n)):
        if any(symp(w,r,n) for r in R):
            continue
        members=[w^r for r in shifts]
        if w!=min(members):
            continue
        cosets+=1
        profile=tuple(sorted(weight(v,n) for v in members))
        if min(profile)>=6:
            deep[str(profile)]+=1
            if profile[0]%2:
                odd[str(profile)]+=1
                assert profile==(7,7,7,7)
    assert cosets==16384
    assert odd=={'(7, 7, 7, 7)':1890}
    K=lambda j,w:sum((-1)**i*3**(j-i)*comb(w,i)*comb(n-w,j-i)
                    for i in range(max(0,j-(n-w)),min(j,w)+1))
    coefficient_checks=0
    # Check both displayed formulas for every nonnegative spectrum of size64.
    for A6 in range(64):
        for A7 in range(64-A6):
            for A8 in range(64-A6-A7):
                A9=63-A6-A7-A8
                n1=K(1,0)+A6*K(1,6)+A7*K(1,7)+A8*K(1,8)+A9*K(1,9)
                n2=K(2,0)+A6*K(2,6)+A7*K(2,7)+A8*K(2,8)+A9*K(2,9)
                assert n1==4*(54-A7-2*A8-3*A9)
                assert n2==16*(A8+3*A9-27)
                if A9==0 and n1>=0 and n2>=0:
                    assert A7==0 and A6==36 and A8==27
                coefficient_checks+=1
    # Explicit size64 even code at the sharp length-nine boundary.
    def multiply(a,b):
        out=0
        while b:
            if b&1:out^=a
            b>>=1;a<<=1
            if a&4:a^=7
        return out
    columns=[(0,1,a) for a in [1,2,3]]+[(1,0,a) for a in [1,2,3]]+[(1,a,0) for a in [1,2,3]]
    generators=[]
    for i in range(3):
        for scalar in [1,2]:
            values=[multiply(scalar,col[i]) for col in columns]
            generators.append(sum((v&1)<<j for j,v in enumerate(values)) |
                              (sum((v>>1)<<j for j,v in enumerate(values))<<9))
    assert len(reduced_basis(generators))==6
    spectrum=Counter(weight(v,9) for v in span(generators))
    assert spectrum=={0:1,6:36,8:27}
    assert len(reduced_basis(gram(generators,9)))==0
    u=next(v for v in span(generators) if weight(v,9)==8)
    values=[multiply(2,((u>>j)&1)|(((u>>(9+j))&1)<<1)) for j in range(9)]
    v=sum((a&1)<<j for j,a in enumerate(values)) | (sum((a>>1)<<j for j,a in enumerate(values))<<9)
    assert v in span(generators) and all(weight(w,9)==8 for w in [u,v,u^v])
    assert ((u&511)|(u>>9))==((v&511)|(v>>9))==(((u^v)&511)|((u^v)>>9))
    files=[Path(__file__),BASE/'model.py',BASE/'RADICAL_PAIR.md']
    report=dict(status='PASS',normalizer_cosets=cosets,deep_profiles=dict(deep),
                odd_deep_profiles=dict(odd),composition_checks=compositions,
                spectrum_coefficient_checks=coefficient_checks,
                sharp_size64_positive_control=dict(spectrum=dict(spectrum),hull=6,
                                                   binary_generators=generators,
                                                   radical_pair=[u,v]),
                sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                scope='Finite checks of a general spectrum restriction and one conditional '
                      'target17 structural exclusion; not a complete parameter result or novelty claim.')
    (BASE/'radical_pair_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))


if __name__=='__main__':
    verify()
