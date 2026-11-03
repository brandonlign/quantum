"""Exact leaf-set filtering plus a physical symmetry/positive control."""
import hashlib
import json
import random
import subprocess
from pathlib import Path
from verify_xdom_prune import run,ExtensionModel
from suffix_symmetries import syndrome_orbits
from signed_extension import solve_extension
from model import span,weight

HERE=Path(__file__).resolve().parent


def restricted(m,d,D,G,reps):
    payload=f'{m} {len(D)} {d}\n'+' '.join(map(str,D))+'\n'+' '.join(map(str,G))+'\n'
    payload+=str(len(reps))+'\n'+' '.join(map(str,reps))+'\n'
    p=subprocess.run([str(HERE/'z_enum_symmetry'),'--gamma','--avoid-radical','--first-orbits'],
                     input=payload,text=True,capture_output=True,check=True,timeout=30)
    stats=json.loads(p.stderr);assert stats['status']=='COMPLETE'
    return [tuple(map(int,line.split()[1:])) for line in p.stdout.splitlines()]


def verify():
    rng=random.Random(411736);cases=[]
    for case in range(80):
        m=3+case%2;b=3+case%3;D=[0]+[rng.randint(1,4) for _ in range((1<<b)-1)];G=[0]*b
        for i in range(b):
            for j in range(i):
                if rng.randrange(2):G[i]^=1<<j;G[j]^=1<<i
        _,L=run('z_enum_avoid_radical',m,3,D,G,['--gamma','--avoid-radical'])
        reps=[v for v in range(1<<b) if D[v]>=2 and rng.randrange(2)]
        actual=restricted(m,3,D,G,reps)
        assert sorted(actual)==sorted(xs for xs in L if xs[0] in reps)
        cases.append(dict(case=case,unrestricted=len(L),restricted=len(actual)))
    model=ExtensionModel(5,0,3,3,[7<<3],[])
    data=dict(labels=list(range(16)),units=[1<<i for i in range(4)],
              distances=model.distances,Gamma=model.Gamma)
    # syndrome_orbits uses threshold five for target17. For this distance-
    # three control obtain all syndrome generators, then build high>=2 orbits.
    orbit=syndrome_orbits(model,data,[]);parent=list(range(16))
    from model import lift_coordinate
    def find(v):
        while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
        return v
    for columns in orbit['verified_syndrome_generators']:
        for v in range(16):
            a,b=find(v),find(lift_coordinate(v,columns))
            if a!=b:parent[max(a,b)]=min(a,b)
    reps=sorted({find(v) for v,d in enumerate(model.distances) if d>=2})
    L=restricted(3,3,model.distances,model.Gamma,reps)
    assert L
    positive=None
    for xs in L:
        a=solve_extension(model,list(xs),minimum_stabilizer_weight=3)
        if a['status']=='FOUND':
            S,D=model.construct(a['model_columns']);assert min(weight(w,5) for w in span(D) if w)==3
            positive=dict(xs=xs,generators=S);break
    assert positive
    out=dict(status='PASS',cases=cases,positive_control=positive,
             physical_automorphism_order=orbit['automorphism_order'],first_representatives=reps,
             sha256={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                     ['z_enum_symmetry.c','z_enum_symmetry','suffix_symmetries.py','verify_first_orbits.py']},
             scope='Complete first-image filtering comparisons and genuine positive code after verified physical automorphisms')
    (HERE/'first_orbit_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print({k:v for k,v in out.items() if k not in ['cases','sha256']})


if __name__=='__main__':verify()
