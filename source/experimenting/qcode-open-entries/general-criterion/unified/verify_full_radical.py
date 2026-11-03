"""Independent finite controls for the full-support radical binary bound."""
from collections import Counter
from math import comb
from pathlib import Path
import hashlib,json


def words(gens):
    out=[0]
    for g in gens:out += [w^g for w in out]
    return out


def wt(w,n):return ((w&((1<<n)-1))|(w>>n)).bit_count()


def sp(u,v,n):
    mask=(1<<n)-1
    return (((u&mask)&(v>>n)).bit_count()+((v&mask)&(u>>n)).bit_count())%2


def main():
    n=9;u=(1<<n)-1
    permutation=[((v<<1)^(11 if v&4 else 0)) for v in range(8)]
    assert sorted(permutation)==list(range(8))
    funcs=[lambda v:1]+[lambda v,i=i:(v>>i)&1 for i in range(3)]
    gens=[u]+[sum(fun(permutation[j])<<j for j in range(8)) |
              (sum(fun(j)<<j for j in range(8))<<n) for fun in funcs]
    C=words(gens);assert len(set(C))==32
    assert all(sp(a,b,n)==0 for a in gens for b in gens)
    spectrum=Counter(wt(w,n) for w in C)
    assert spectrum=={0:1,6:14,7:14,8:1,9:2}
    image=sorted(set(w>>n for w in C));assert len(image)==16
    assert all(v.bit_count()%2==0 for v in image)
    assert min(v.bit_count() for v in image if v)==4
    assert {w for w in C if w>>n==0}=={0,u}
    for w in C:
        if w not in [0,u]:assert wt(w,n)+wt(w^u,n)==n+(w>>n).bit_count()
    # Exactly the admissible rank multiplicities and zero allocations.
    geometries=[]
    for rank0 in range(2):
        for z1 in range(2):
            for z2 in range(2):
                for z3 in range(2):
                    if any(rank0+z>1 for z in [z1,z2,z3]):continue
                    rank1=z1+z2+z3;rank2=9-rank0-rank1
                    if rank2<0 or rank2%2:continue
                    geometries.append((rank0,rank1,rank2))
    assert set(geometries)=={(1,0,8),(0,1,8),(0,3,6)}
    pairs={'common_zero':(255,255<<9),'one_zero':(511,255<<9),
           'different_zeros':(510,5|(504<<9))}
    deep={}
    for name,(a,b) in pairs.items():
        assert sp(a,b,9)==0
        assert min(wt(v,9) for v in [a,b,a^b])>=8
        profiles=Counter();cosets=0
        for w in range(1<<18):
            if sp(w,a,9) or sp(w,b,9):continue
            orbit=[w,w^a,w^b,w^a^b]
            if w!=min(orbit):continue
            cosets+=1;profile=tuple(sorted(wt(v,9) for v in orbit))
            if min(profile)>=6:profiles[str(profile)]+=1
        assert cosets==16384
        deep[name]=dict(generators=[a,b],weights=[wt(v,9) for v in [a,b,a^b]],
                        normalizer_cosets=cosets,deep_profiles=dict(profiles))
    assert (1<<5)*sum(comb(8,i) for i in range(2))==288>256
    report=dict(status='PASS',sharp_control=dict(binary_generators=gens,
        spectrum=dict(spectrum),hull=5,binary_projection_dimension=4,
        binary_projection_distance=4,full_support_radical=u),
        rank_geometries=sorted(set(geometries)),canonical_suffix_pairs=deep,
        scope='Conditional additive bound and exact heavy suffix-stabilizer geometry reduction; '
              'not a complete target17 exclusion or originality claim.',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [Path(__file__),Path(__file__).with_name('FULL_RADICAL_BOUND.md')]})
    Path(__file__).with_name('full_radical_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
