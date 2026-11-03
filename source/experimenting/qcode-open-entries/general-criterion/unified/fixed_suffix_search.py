"""Exact fixed-suffix quotient input for rank-one localized extension searches.

The C enumerator has generic m,d input. Exact-N pruning requires independent
Z images: R=<Z^m>, d>floor(m/2), and shortened-code minimum distance>=d
(including zero). Low-weight stabilizers confined to R are permitted.
This module does not assume fixing one N covers a whole class.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from model import ExtensionModel,coordinates,gram,lift_coordinate,reduced_basis,span,weight,nullspace
from run_benchmark import z_condition
from signed_extension import solve_extension

HERE=Path(__file__).resolve().parent


def quotient_input(model,N):
    N=reduced_basis(N);ell=len(N)
    assert model.r==1 and model.R==[((1<<model.m)-1)<<model.m]
    assert model.d>model.m//2 and (model.C_distance is None or model.C_distance>=model.d)
    Ncoords=[coordinates(v,model.P) for v in N]
    assert all(gram(model.P+[v],model.f)[-1]==0 for v in N)
    full=list(Ncoords);quotient=[]
    for i in range(model.b):
        v=1<<i
        if len(reduced_basis(full+[v]))>len(full):full.append(v);quotient.append(v)
    assert len(full)==model.b
    units=[sum(((coordinates(1<<i,full)>>(ell+j))&1)<<i for i in range(model.b))
           for j in range(model.b-ell)]
    labels=[lift_coordinate(v,units) for v in range(1<<(model.b-ell))]
    assert len(set(labels))==len(labels)
    assert all(all((old&n).bit_count()%2==0 for n in Ncoords) for old in labels)
    Q=[lift_coordinate(v,model.P) for v in quotient]
    return dict(units=units,labels=labels,distances=[model.distances[v] for v in labels],Gamma=gram(Q,model.f))


def search(model,N,minimum_stabilizer_weight,seconds=30,solver_nodes=10000):
    data=quotient_input(model,N);binary=HERE/'m8r1/z_enum_avoid_radical'
    payload=f'{model.m} {len(data["labels"])} {model.d}\n'+' '.join(map(str,data['distances']))+'\n'
    payload+=' '.join(map(str,data['Gamma']))+'\n'
    tally=Counter();found=None;leaves=0
    with tempfile.TemporaryFile('w+') as err:
        proc=subprocess.Popen([str(binary),'--gamma','--avoid-radical','--seconds',str(seconds)],
                              stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True)
        try:
            proc.stdin.write(payload);proc.stdin.close()
            for line in proc.stdout:
                fields=line.split();assert fields[0]=='LEAF' and len(fields)==model.m
                xs=[lift_coordinate(int(v),data['units']) for v in fields[1:]];leaves+=1
                assert z_condition(xs,model.Gamma,model.b)
                answer=solve_extension(model,xs,node_limit=solver_nodes,
                                       minimum_stabilizer_weight=minimum_stabilizer_weight)
                tally[answer['status']+':'+answer.get('reason','')]+=1
                if answer['status']=='FOUND':found=dict(xs=xs,answer=answer);break
            if found:proc.terminate()
            proc.wait();err.seek(0)
            if found:stats=dict(status='FOUND_EARLY',leaves=leaves)
            else:
                stats=json.loads(err.read());assert proc.returncode==0
                assert stats['status'] in ['COMPLETE','LIMIT'] and stats['leaves']==leaves
        finally:
            if proc.poll() is None:proc.terminate();proc.wait()
    status='FOUND' if found else ('EXCLUDED_FIXED_N' if stats['status']=='COMPLETE' and
             not any(v for k,v in tally.items() if not k.startswith('EXCLUDED:')) else 'LIMIT')
    return dict(status=status,z=stats,tally=dict(tally),found=found,N=list(N),
                quotient_dimension=model.b-len(reduced_basis(N)),
                scope='Only extensions with exactly this suffix stabilizer; class exclusion needs complete justified N coverage')


def main():
    witness_path=HERE/'suffix10_radical_pair_witness.json'
    witness=json.loads(witness_path.read_text());C=witness['generators']
    model=ExtensionModel(17,3,6,7,[127<<7],C)
    assert model.hull==3
    radical=[lift_coordinate(v,C) for v in nullspace(gram(C,10),len(C))]
    heavy={w for w in span(radical) if weight(w,10)>=7}
    planes=[(a,b) for a in sorted(heavy) for b in sorted(heavy) if a<b<(a^b) and a^b in heavy]
    assert len(planes)==1
    result=search(model,planes[0],7)
    result.update(parameters=dict(n=17,k=3,d=6,m=7,r=1),heavy_plane_count=1,
                  sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                          [witness_path,Path(__file__),HERE/'model.py',HERE/'signed_extension.py',
                           HERE/'m8r1/z_enum_avoid_radical.c',HERE/'m8r1/z_enum_avoid_radical']})
    (HERE/'suffix10_witness_extension.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['sha256','found']})


if __name__=='__main__':main()
