"""All minimum-six N planes in the recentered rank-one hull-five cell."""
import json
from pathlib import Path
import sys
from collections import deque

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import lift_coordinate, reduced_basis, weight, symp
from suffix_symmetries import automorphisms
from mixed_frames.stage3.common import atomic_json,sha

SOURCE=ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'


def plane_key(a,b):return tuple(sorted((a,b,a^b)))


def prepare():
    classes=json.loads(SOURCE.read_text())['classes']
    records=[]
    for c in classes:
        if c['representative']['hull']!=5:continue
        C=reduced_basis(c['representative']['generators'])
        assert len(C)==5
        words=[lift_coordinate(i,C) for i in range(32)]
        assert min(weight(w,9) for w in words if w)>=6
        assert all(symp(a,b,9)==0 for a in C for b in C)
        planes={plane_key(words[a],words[b]) for a in range(1,32) for b in range(a+1,32)}
        assert len(planes)==155
        maps,order=automorphisms(C,9)
        remaining=set(planes);orbits=[]
        while remaining:
            seed=min(remaining);seen={seed};todo=deque([seed])
            while todo:
                p=todo.popleft()
                for apply in maps:
                    q=plane_key(apply(p[0]),apply(p[1]))
                    assert q in planes
                    if q not in seen:seen.add(q);todo.append(q)
            assert seen<=remaining
            remaining-=seen
            orbits.append(dict(representative=list(seed[:2]),members=[list(p) for p in sorted(seen)]))
        assert sum(len(o['members']) for o in orbits)==155
        records.append(dict(cls=c['class_index'],C=C,plane_count=155,
                            automorphism_order=order,orbits=orbits))
    assert len(records)==118
    result=dict(status='PASS',minimum_stabilizer_weight=6,records=records,
                classes=118,physical_planes=18290,plane_orbits=sum(len(r['orbits']) for r in records),
                source_sha256=sha(SOURCE),preparation_sha256=sha(Path(__file__)),
                scope='Complete N-plane partition of the classified hull-five cell under minimum six; not an exclusion')
    atomic_json(HERE/'plane_plan.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)


if __name__=='__main__':prepare()
