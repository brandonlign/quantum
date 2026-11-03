"""Physical quotient checks, a real-code positive control, and alternate m7 replay."""
import hashlib
import json
from pathlib import Path
import subprocess
from model import ExtensionModel,span,weight,symp,split
from fixed_suffix_search import quotient_input,search

HERE=Path(__file__).resolve().parent


def main():
    path=HERE/'radical_pair_lift10.json';report=json.loads(path.read_text())
    assert report['status']=='COMPLETE_FIXED_N_SUBFAMILY_EXCLUSION'
    for name,digest in report['sha256'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest
    records=[];cosets=0
    for r in report['records']:
        model=ExtensionModel(17,3,6,7,[127<<7],r['generators']);q=quotient_input(model,r['N'])
        assert model.hull==3 and model.C_distance==6
        # Independent physical minimization for every restricted coset.
        words=span(model.C)
        for v,old in enumerate(q['labels']):
            w=model.representatives[old]
            assert q['distances'][v]==min(weight(w^c,10) for c in words)
            cosets+=1
        payload=f'7 {len(q["labels"])} 6\n'+' '.join(map(str,q['distances']))+'\n'+' '.join(map(str,q['Gamma']))+'\n'
        # This backend computes full annihilator Gram matrices and omits
        # the exact-N radical-disjointness prune used in the primary run.
        run=subprocess.run([str(HERE/'m8r1/z_enum_xdom'),'--gamma','--seconds','10'],
                           input=payload,text=True,capture_output=True,check=True)
        stats=json.loads(run.stderr)
        assert stats['status']=='COMPLETE' and stats['leaves']==0 and not run.stdout.strip()
        records.append(dict(orbit=r['orbit'],z=stats))
    control=ExtensionModel(7,1,3,4,[15<<4],[7]);positive=search(control,[],3,seconds=10)
    assert positive['status']=='FOUND'
    columns=positive['found']['answer']['model_columns'];assert control.criterion(columns)
    S,D=control.construct(columns)
    assert all(symp(a,b,7)==0 for a in S for b in S)
    assert min(weight(w,7) for w in span(D) if w)==3
    stabilizer_minimum=min(weight(w,7) for w in span(S) if w)
    assert stabilizer_minimum>=3
    assert not any(w for w in span(S) if split(w,4,3)[0]==0)
    zero=ExtensionModel(5,0,3,3,[7<<3],[])
    zero_result=search(zero,[],3,seconds=10)
    assert zero_result['status']=='FOUND'
    zero_columns=zero_result['found']['answer']['model_columns']
    assert zero.criterion(zero_columns)
    zeroS,zeroD=zero.construct(zero_columns)
    assert all(symp(a,b,5)==0 for a in zeroS for b in zeroS)
    assert min(weight(w,5) for w in span(zeroD) if w)==3
    assert not any(w for w in span(zeroS) if split(w,3,2)[0]==0)
    result=dict(status='PASS',restricted_physical_cosets=cosets,alternate_complete_orbits=len(records),
                alternate_nodes=sum(r['z']['nodes'] for r in records),records=records,
                positive_control=dict(parameters=[7,1,3],prefix_length=4,generators=S,
                                      pure_distance=3,minimum_stabilizer_weight=stabilizer_minimum,
                                      suffix_stabilizer_dimension=0),
                zero_suffix_positive_control=dict(parameters=[5,0,3],prefix_length=3,
                                                  generators=zeroS,pure_distance=3),
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [path,Path(__file__),HERE/'fixed_suffix_search.py',HERE/'model.py',
                         HERE/'signed_extension.py',HERE/'m8r1/z_enum_xdom.c',HERE/'m8r1/z_enum_xdom']},
                scope='All restricted cosets physically checked and ten complete alternate Gram-prune replays; canonical distance core shared; no second complete lift classification')
    (HERE/'fixed_suffix_search_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['records','sha256','positive_control']})


if __name__=='__main__':main()
