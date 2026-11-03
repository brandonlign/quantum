"""Independent physical checks of the label/subcode inverse characterization."""
from pathlib import Path
import json,random,hashlib
from gauge_fix import fix_prefix_kernel
from subcode_selection import SelectionContext
from model import reduced_basis,span,symp,weight
BASE=Path(__file__).resolve().parent


def main():
    rng=random.Random(175514);checks=0;distance_assertions=0;minimum_assertions=0;lowered=0;first_counterexample=None
    for n in range(2,6):
        for trial in range(20):
            rank=rng.randrange(n+1);S=[1<<(n+i) for i in range(rank)]
            for step in range(20):
                i=rng.randrange(n);j=rng.choice([v for v in range(n) if v!=i]);typ=rng.randrange(3)
                if typ==0:S=[v ^ ((((v>>i)&1)^((v>>(n+i))&1))*((1<<i)|(1<<(n+i)))) for v in S]
                elif typ==1:S=[v^(((v>>i)&1)<<(n+i)) for v in S]
                else:S=[v^(((v>>i)&1)<<j)^(((v>>(n+j))&1)<<(n+i)) for v in S]
            T=rng.sample(range(n),rng.randrange(n+1));fixed=fix_prefix_kernel(S,n,T);R=fixed['R'];ell=fixed['fixed_logicals']
            oldbasis=list(R);complement=[]
            for v in reduced_basis(S):
                if len(reduced_basis(oldbasis+[v]))>len(oldbasis):oldbasis.append(v);complement.append(v)
            U=list(R)+complement+fixed['L'];context=SelectionContext(U,n)
            assert len(U)==len(fixed['S_fixed'])
            for choice in range(3):
                rows=[]
                for j in range(ell):
                    row=1<<(len(oldbasis)+j)
                    for i in range(len(R),len(oldbasis)):
                        if rng.randrange(2):row |= 1<<i
                    rows.append(row)
                # Direct Pauli membership, not syndrome shortest-path logic.
                evaluated=context.evaluate(rows,1,1,T,R);Sn=evaluated['generators'];st=set(span(Sn))
                D={w for w in range(1<<(2*n)) if all(symp(w,s,n)==0 for s in Sn)}
                actual_d=min((weight(w,n) for w in D-st),default=None)
                actual_minS=min((weight(w,n) for w in st if w),default=None)
                for d in range(1,n+2):
                    r=context.evaluate(rows,d,1,T,R)
                    assert r['distance_ok']==(actual_d is None or actual_d>=d)
                    assert r['local_and_shortened_code_ok'];distance_assertions+=1
                for mu in range(1,n+2):
                    r=context.evaluate(rows,1,mu,T,R)
                    assert r['minimum_stabilizer_ok']==(actual_minS is None or actual_minS>=mu);minimum_assertions+=1
                if actual_d is not None and (context.base_distance is None or context.base_distance>actual_d):
                    lowered+=1
                    if first_counterexample is None:first_counterexample=dict(n=n,base_generators=U,label_rows=rows,
                        base_distance=context.base_distance,subcode_distance=actual_d,
                        syndrome_rowcode_distance=evaluated['nonzero_syndrome_rowcode_distance'])
                checks+=1
    # High-distance positive base: fixing one logical of the five-qubit code.
    pauli=lambda text:sum((c in 'XY')<<i for i,c in enumerate(text)) | (sum((c in 'YZ')<<i for i,c in enumerate(text))<<len(text))
    S5=[pauli(t) for t in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']]
    fixed=fix_prefix_kernel(S5,5,[0,1,2,3]);R=fixed['R'];old=list(R)
    for v in reduced_basis(S5):
        if len(reduced_basis(old+[v]))>len(old):old.append(v)
    U=old+fixed['L'];context=SelectionContext(U,5);five=[]
    for bits in range(1<<(len(old)-len(R))):
        row=(1<<len(old))|(bits<<len(R));r=context.evaluate([row],3,3,[0,1,2,3],R)
        five.append(dict(label_row=row,distance_three=r['distance_ok'],minimum_stabilizer_three=r['minimum_stabilizer_ok'],
                         syndrome_coset_distance=r['nonzero_syndrome_rowcode_distance']))
    assert any(v['distance_three'] for v in five) and any(not v['distance_three'] for v in five)
    report=dict(status='PASS',subcode_selections=checks,distance_conditions_checked=distance_assertions,
        minimum_stabilizer_conditions_checked=minimum_assertions,selections_lowering_distance=lowered,
        inverse_distance_counterexample=first_counterexample,five_qubit_high_distance_base=dict(generators=U,minimum_stabilizer_weight=min(context.weights[1:]),selections=five),
        scope='Independent physical-membership distances and stabilizer weights for all tested labels; '
              'the theorem is general, this implementation is a bounded reference.',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [Path(__file__),BASE/'subcode_selection.py',BASE/'gauge_fix.py',BASE/'model.py']})
    (BASE/'subcode_selection_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
