"""Verified suffix automorphism orbits of heavy stabilizer planes."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from run_hull5_planes import heavy_planes,reduced_basis,weight
from suffix_symmetries import automorphisms

HERE=Path(__file__).resolve().parent


def key(a,b):return tuple(sorted([a,b,a^b]))


def geometry(a,b):
    words=[a,b,a^b]
    if any(weight(w,9)==9 for w in words):return 'full_support_word'
    support=[(w&511)|(w>>9) for w in words]
    if (~support[0])&(~support[1])&(~support[2])&511:return 'common_zero'
    return 'different_zeros'


def main():
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    classes=json.loads(source.read_text())['classes'];records=[]
    for c in classes:
        if c['representative']['hull']!=5:continue
        C=reduced_basis(c['representative']['generators']);planes=heavy_planes(C)
        keys=[key(a,b) for a,b in planes];lookup={k:i for i,k in enumerate(keys)}
        parent=list(range(len(planes)))
        def find(v):
            while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
            return v
        maps,order=automorphisms(C,9)
        for apply in maps:
            for i,(a,b) in enumerate(planes):
                j=lookup[key(apply(a),apply(b))]
                assert geometry(a,b)==geometry(*planes[j])
                x,y=find(i),find(j)
                if x!=y:parent[max(x,y)]=min(x,y)
        groups={}
        for i in range(len(planes)):groups.setdefault(find(i),[]).append(i)
        orbits=[dict(representative=list(planes[min(indices)]),geometry=geometry(*planes[min(indices)]),
                     plane_indices=indices,size=len(indices)) for indices in groups.values()]
        records.append(dict(cls=c['class_index'],automorphism_order=order,plane_count=len(planes),orbits=orbits))
    result=dict(status='PASS',records=records,physical_planes=sum(r['plane_count'] for r in records),
                orbit_count=sum(len(r['orbits']) for r in records),
                remaining_different_zero_orbits=sum(o['geometry']=='different_zeros' for r in records for o in r['orbits']),
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [source,HERE/'suffix_symmetries.py',Path(__file__)]},
                scope='Automorphisms physically preserve C, Pauli weights, symplectic form and plane geometry; one representative per orbit suffices for extension existence')
    (HERE/'heavy_plane_orbits.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['records','sha256']})


if __name__=='__main__':main()
