"""Physical-membership and positive-control checks for prefix kernel fixing."""
from pathlib import Path
from collections import Counter
import json,random,hashlib
from gauge_fix import fix_prefix_kernel
from model import span,weight,reduced_basis,symp
BASE=Path(__file__).resolve().parent


def pauli(text):
    n=len(text)
    return sum((c in 'XY')<<i for i,c in enumerate(text)) | (sum((c in 'YZ')<<i for i,c in enumerate(text))<<n)


def main():
    rng=random.Random(1601755);checks=0;types=Counter()
    for n in range(2,7):
        for _ in range(40):
            r=rng.randrange(n+1);S=[1<<(n+i) for i in range(r)]
            for gate in range(30):
                typ=rng.randrange(3);i=rng.randrange(n)
                if typ==0:
                    S=[v ^ ((((v>>i)&1)^((v>>(n+i))&1))*((1<<i)|(1<<(n+i)))) for v in S]
                elif typ==1:S=[v^(((v>>i)&1)<<(n+i)) for v in S]
                else:
                    j=rng.choice([a for a in range(n) if a!=i]);S=[v^(((v>>i)&1)<<j)^(((v>>(n+j))&1)<<(n+i)) for v in S]
            T=rng.sample(range(n),rng.randrange(n+1));res=fix_prefix_kernel(S,n,T)
            St=set(span(S));Sf=set(span(res['S_fixed']));D={w for w in range(1<<(2*n)) if all(symp(w,s,n)==0 for s in S)}
            Df={w for w in D if all(symp(w,s,n)==0 for s in res['S_fixed'])}
            Tmask=sum((1<<i)|(1<<(n+i)) for i in T);full=(1<<(2*n))-1
            R={w for w in St if not w&(full^Tmask)};Rf={w for w in Sf if not w&(full^Tmask)}
            C={w for w in D if not w&Tmask};Cf={w for w in Df if not w&Tmask}
            assert C==Cf==set(span(res['shortened_centralizer']))
            assert Rf=={w for w in Df if not w&(full^Tmask)}
            assert len(Sf)==len(St)*(1<<res['fixed_logicals'])
            if Df-Sf:assert min(weight(w,n) for w in Df-Sf)>=min(weight(w,n) for w in D-St)
            checks+=1;types[res['kernel_dimension'],res['kernel_form_rank']]+=1
    S5=[pauli(s) for s in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']];five=fix_prefix_kernel(S5,5,[0,1,2,3])
    root=BASE.parents[1];S16=[]
    for line in (root/'release16/search/q1645.txt').read_text().splitlines():
        x,z=line.split('|');S16.append(sum(int(v)<<i for i,v in enumerate(x.split())) | (sum(int(v)<<i for i,v in enumerate(z.split()))<<16))
    p=min((w for w in span(S16) if w),key=lambda w:weight(w,16));T=[i for i in range(16) if ((p&65535)|(p>>16))>>i&1]
    sixteen=fix_prefix_kernel(S16,16,T)
    compact=lambda v:{k:v[k] for k in ['n','k','k_fixed','local_rank','kernel_dimension','kernel_form_rank','fixed_logicals','local_rank_fixed','fixed_kernel_dimension']}
    report=dict(status='PASS',small_physical_checks=checks,kernel_types=[dict(q=q,t=t,count=c) for (q,t),c in sorted(types.items())],
        five_qubit_control=compact(five),sixteen_qubit_control=compact(sixteen),
        scope='Direct physical membership checks for small codes, and algebraic checks on known larger controls; '
              'no full large-code distance rerun or novelty claim.',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),BASE/'gauge_fix.py',BASE/'model.py']})
    (BASE/'gauge_fix_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
