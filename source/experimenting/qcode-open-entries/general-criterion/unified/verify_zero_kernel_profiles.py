"""Direct physical validation of fixed radical distance tables, including impurity."""
from pathlib import Path
import random,json,hashlib
from gauge_fix import fix_prefix_kernel
from zero_kernel_profiles import recover_profiles
from model import span,weight,symp,reduced_basis,embed
BASE=Path(__file__).resolve().parent


def direct(S,n):
    st=set(span(S));D={w for w in range(1<<(2*n)) if all(symp(w,s,n)==0 for s in S)}
    return dict(logical_distance=min((weight(w,n) for w in D-st),default=None),
                minimum_stabilizer_weight=min((weight(w,n) for w in st if w),default=None),
                pure_distance=min((weight(w,n) for w in D if w),default=None))


def main():
    rng=random.Random(175506);checks=0;states=0;impure=0
    for n in range(2,7):
        for _ in range(20):
            r=rng.randrange(n+1);S=[1<<(n+i) for i in range(r)]
            for step in range(20):
                i=rng.randrange(n);j=rng.choice([v for v in range(n) if v!=i]);typ=rng.randrange(3)
                if typ==0:S=[v ^ ((((v>>i)&1)^((v>>(n+i))&1))*((1<<i)|(1<<(n+i)))) for v in S]
                elif typ==1:S=[v^(((v>>i)&1)<<(n+i)) for v in S]
                else:S=[v^(((v>>i)&1)<<j)^(((v>>(n+j))&1)<<(n+i)) for v in S]
            T=rng.sample(range(n),rng.randrange(n+1));fixed=fix_prefix_kernel(S,n,T);Sf=fixed['S_fixed']
            model,columns=recover_profiles(Sf,n,T);actual=direct(Sf,n)
            assert model.metrics(columns)==actual
            states+=actual['logical_distance'] is None
            impure+=actual['logical_distance'] is not None and actual['minimum_stabilizer_weight'] is not None and actual['minimum_stabilizer_weight']<actual['logical_distance']
            checks+=1
    pauli=lambda text:sum((c in 'XY')<<i for i,c in enumerate(text)) | (sum((c in 'YZ')<<i for i,c in enumerate(text))<<len(text))
    S5=[pauli(t) for t in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']]
    S6=[embed(pauli('Z'),0,1,5)]+[embed(0,w,1,5) for w in S5]
    S7=[embed(pauli('Z'),0,1,6)]+[embed(0,w,1,6) for w in S6]
    controls=[]
    for S,n in [(S6,6),(S7,7)]:
        model,columns=recover_profiles(S,n,[0]);metrics=model.metrics(columns);assert metrics==direct(S,n)
        assert metrics['logical_distance']==3 and metrics['minimum_stabilizer_weight']==1
        controls.append(dict(n=n,k=n-len(S),suffix_code_distance=model.min_C,profiles=metrics))
    assert controls[1]['suffix_code_distance']==1
    result=dict(status='PASS',small_physical_codes=checks,state_controls=states,impure_random_controls=impure,
        explicit_impure_controls=controls,
        scope='Direct Pauli-membership logical distances and stabilizer/pure weights; '
              'fixed relative profiles accept a real nonlocalized impure case with C distance1.',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [Path(__file__),BASE/'zero_kernel_profiles.py',BASE/'gauge_fix.py',BASE/'model.py']})
    (BASE/'zero_kernel_profile_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='sha256'}))
if __name__=='__main__':main()
