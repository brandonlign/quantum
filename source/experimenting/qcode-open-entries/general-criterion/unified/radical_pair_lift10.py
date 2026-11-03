"""Complete two-site radical-plane lifts with an eight-site rank-five kernel.

Enumerates a structural subfamily, not the whole (10,128,6) cell.
The equivalence argument also allows simultaneous S3 on the two tail sites.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'m8r1'))
from suffix_symmetries import automorphisms
from model import span,dual,reduced_basis,weight,gram,embed,coordinates,ExtensionModel
from fixed_suffix_search import search

HERE=Path(__file__).resolve().parent


def main():
    source=HERE/'critical_prefix_spectrum_verification.json'
    G=json.loads(source.read_text())['positive_control_generators'];words=span(G)
    fullwords=[w for w in words if weight(w,8)==8]
    assert len(fullwords)==3 and fullwords[0]^fullwords[1]==fullwords[2]
    N=fullwords[:2];B=reduced_basis(G);Q=[]
    for v in dual(N,8):
        if len(reduced_basis(B+[v]))>len(B):Q.append(v);B=reduced_basis(B+[v])
    Q=reduced_basis(Q);full=G+Q;reps=span(Q)
    assert len(Q)==9 and all(coordinates(w,full)>>5==i for i,w in enumerate(reps))
    delta=[min(weight(v^w,8) for w in words) for v in reps]
    assert max(delta)==4
    heavy=[i for i,d in enumerate(delta) if d>=4]
    def code(a,b):
        return [embed(v,0,8,2) for v in G]+[embed(reps[a],3,8,2),embed(reps[b],12,8,2)]
    keys=set();labelled_pairs=0
    for a in heavy:
        for b in heavy:
            if b<=a or delta[a^b]<4:continue
            if len(reduced_basis(gram(code(a,b),10)))==4:
                keys.add(tuple(sorted([a,b,a^b])));labelled_pairs+=1
    assert labelled_pairs==3696 and len(keys)==1232
    maps,order=automorphisms(G,8)
    transforms=[[coordinates(g(v),full)>>5 for v in reps] for g in maps]
    assert all(sorted(T)==list(range(512)) and all(delta[i]==delta[T[i]] for i in range(512)) for T in transforms)
    remaining=set(keys);orbits=[]
    while remaining:
        seed=min(remaining);orbit={seed};front=[seed]
        while front:
            pair=front.pop()
            for T in transforms:
                image=tuple(sorted(T[v] for v in pair));assert image in keys
                if image not in orbit:orbit.add(image);front.append(image)
        assert orbit<=remaining;remaining-=orbit
        orbits.append(dict(seed=list(seed),size=len(orbit)))
    assert len(orbits)==10 and sum(o['size'] for o in orbits)==1232
    fixedN=[embed(v,0,8,2) for v in N];records=[]
    for i,o in enumerate(orbits):
        C=code(*o['seed'][:2]);assert len(reduced_basis(C))==7
        assert min(weight(w,10) for w in span(C) if w)==6
        assert all(v in span(dual(C,10)) for v in fixedN)
        model=ExtensionModel(17,3,6,7,[127<<7],C);assert model.hull==3
        r=search(model,fixedN,7);r.update(orbit=i,orbit_size=o['size'],generators=C,seed=o['seed'])
        records.append(r);print(json.dumps(dict(orbit=i,status=r['status'],z=r['z'])),flush=True)
    complete=all(r['status']=='EXCLUDED_FIXED_N' for r in records)
    result=dict(status='COMPLETE_FIXED_N_SUBFAMILY_EXCLUSION' if complete else 'PARTIAL_SUBFAMILY',
                labelled_pairs=labelled_pairs,quotient_planes=len(keys),orbit_count=len(orbits),
                prefix_automorphism_order=order,records=records,
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [source,Path(__file__),HERE/'fixed_suffix_search.py',HERE/'model.py',
                         HERE/'signed_extension.py',HERE/'m8r1/suffix_symmetries.py',
                         HERE/'m8r1/z_enum_avoid_radical.c',HERE/'m8r1/z_enum_avoid_radical']},
                scope='Only m7,r1,hull3 with prescribed common-zero eight-site suffix stabilizer and rank-five kernel after shortening its two zero sites; not the whole rank-one branch or target')
    (HERE/'radical_pair_lift10.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['records','sha256']},flush=True)


if __name__=='__main__':main()
