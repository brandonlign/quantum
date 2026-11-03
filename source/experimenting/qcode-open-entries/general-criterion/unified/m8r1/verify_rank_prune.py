"""Complete comparisons for the optimized rank-budget backend."""
import hashlib
import json
import random
from pathlib import Path
from verify_xdom_prune import run, unary_ok, z_condition, ExtensionModel

HERE=Path(__file__).resolve().parent
rng=random.Random(2601736)
cases=[]
for case in range(100):
    m=3+case%2; b=3+case%3
    D=[0]+[rng.randint(1,4) for _ in range((1<<b)-1)]
    G=[0]*b
    for i in range(b):
        for j in range(i):
            if rng.randrange(2):G[i]^=1<<j;G[j]^=1<<i
    s0,L0=run('z_enum8192',m,3,D,G,[])
    expected=sorted(xs for xs in L0 if z_condition(xs,G,b) and unary_ok(xs,D,m,3))
    s1,L1=run('z_enum_rank',m,3,D,G,['--gamma','--xdom'])
    assert sorted(L1)==expected,(case,expected,L1)
    s2,L2=run('z_enum_xdom',m,3,D,G,['--gamma','--xdom'])
    assert sorted(L2)==sorted(L1)
    cases.append(dict(case=case,plain=len(L0),survivors=len(L1),
                      old_nodes=s2['nodes'],new_nodes=s1['nodes']))

# Exercise 13-bit labels without an enormous test tree. Extend a physical
# positive profile by zero-distance labels which cannot occur as Z images.
model=ExtensionModel(5,0,3,3,[7<<3],[])
D=list(model.distances)+[0]*(8192-16)
G=list(model.Gamma)+[0]*9
s0,L0=run('z_enum8192',3,3,D,G,[])
s1,L1=run('z_enum_rank',3,3,D,G,['--gamma','--xdom'])
expected=sorted(xs for xs in L0 if z_condition(xs,G,13) and unary_ok(xs,D,3,3))
assert L1 and sorted(L1)==expected
assert (3,12) in L1
result=dict(status='PASS',instances=cases,
            nonempty_cases=sum(c['survivors']>0 for c in cases),
            label8192_control=dict(leaves=len(L1),retains_physical_positive=True),
            sha256={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                    for n in ['z_enum_rank.c','z_enum_rank','verify_rank_prune.py']},
            scope='Exact complete leaf sets, independent physical-weight oracle, old backend, and 8192-label control')
(HERE/'rank_prune_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['instances','sha256']}))
