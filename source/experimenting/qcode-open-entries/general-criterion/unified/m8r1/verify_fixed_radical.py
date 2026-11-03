"""Independent final-radical and exhaustive-leaf checks for early pruning."""
import hashlib
import json
import random
from pathlib import Path
from verify_xdom_prune import run,z_condition,ExtensionModel
from model import nullspace,reduced_basis
from model import span,weight,symp
from signed_extension import solve_extension

HERE=Path(__file__).resolve().parent
rng=random.Random(311736);cases=[]
for case in range(120):
    m=3+case%2;b=3+case%3
    D=[0]+[rng.randint(1,4) for _ in range((1<<b)-1)]
    G=[0]*b
    for i in range(b):
        for j in range(i):
            if rng.randrange(2):G[i]^=1<<j;G[j]^=1<<i
    h=b-len(reduced_basis(G))
    _,plain=run('z_enum8192',m,3,D,G,[])
    expected=sorted(xs for xs in plain if z_condition(xs,G,b)
                    and not nullspace(list(G)+list(xs),b))
    _,primary=run('z_enum_fixed_radical',m,3,D,G,['--gamma','--radical-rank',str(h)])
    _,alternate=run('z_enum_avoid_radical',m,3,D,G,['--gamma','--avoid-radical'])
    assert sorted(primary)==expected,(case,'primary',expected,primary)
    assert sorted(alternate)==expected,(case,'alternate',expected,alternate)
    cases.append(dict(case=case,plain=len(plain),survivors=len(primary),radical_dimension=h))
model=ExtensionModel(5,0,3,3,[7<<3],[])
_,primary=run('z_enum_fixed_radical',3,3,model.distances,model.Gamma,['--gamma','--radical-rank','0'])
_,alternate=run('z_enum_avoid_radical',3,3,model.distances,model.Gamma,['--gamma','--avoid-radical'])
assert (3,12) in primary and sorted(primary)==sorted(alternate)
model7=ExtensionModel(7,1,3,4,[15<<4],[7])
_,primary7=run('z_enum_fixed_radical',4,3,model7.distances,model7.Gamma,['--gamma','--radical-rank','1'])
_,alternate7=run('z_enum_avoid_radical',4,3,model7.distances,model7.Gamma,['--gamma','--avoid-radical'])
assert sorted(primary7)==sorted(alternate7) and (3,5,10) in primary7
answer=solve_extension(model7,[3,5,10],minimum_stabilizer_weight=3)
assert answer['status']=='FOUND'
S7,D7=model7.construct(answer['model_columns'])
assert len(S7)==6 and all(symp(a,b,7)==0 for a in S7 for b in S7)
assert min(weight(w,7) for w in span(D7) if w)==3
assert not nullspace(model7.Gamma+[3,5,10],5)
out=dict(status='PASS',cases=cases,nonempty_cases=sum(c['survivors']>0 for c in cases),
         positive_code_retained=True,
         nonzero_radical_positive_control=dict(parameters=[7,1,3],Gamma_radical_dimension=1,
             Z_columns=[3,5,10],generators=S7,pure_distance=3,forced_suffix_stabilizer_dimension=0),
         sha256={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                 ['z_enum_fixed_radical.c','z_enum_fixed_radical','z_enum_avoid_radical.c','z_enum_avoid_radical','verify_fixed_radical.py']},
         scope='Complete oracle comparisons; rank-budget implementation versus full-annihilator/radical implementation; distance enumeration core shared')
(HERE/'fixed_radical_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k not in ['cases','sha256']})
