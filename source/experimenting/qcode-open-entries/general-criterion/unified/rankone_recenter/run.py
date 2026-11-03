"""Resumable bounded rank-one minimum-six screening. LIMIT never excludes."""
import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel, lift_coordinate, span, symp, weight
from run_benchmark import z_condition
from signed_extension import solve_extension
from run_hull5_planes import quotient_input
from fixed_suffix_search import quotient_input as generic_quotient
from suffix_symmetries import syndrome_orbits
from mixed_frames.stage3.common import atomic_json,sha
from prepare import SOURCE

BINARY=ROOT/'m8r1/z_enum_symmetry'


def fingerprints():
    paths=[SOURCE,Path(__file__),HERE/'prepare.py',HERE/'plane_plan.json',BINARY,
           ROOT/'m8r1/z_enum_symmetry.c',ROOT/'m8r1/suffix_symmetries.py',
           ROOT/'m8r1/run_hull5_planes.py',ROOT/'fixed_suffix_search.py',ROOT/'model.py',
           ROOT/'signed_extension.py',ROOT/'kernel_bound.py',ROOT/'two_site.py',ROOT/'run_benchmark.py']
    return {str(p.relative_to(ROOT)):sha(p) for p in paths}


@lru_cache(None)
def target(C):return ExtensionModel(17,3,6,8,[255<<8],list(C))


def enumerate_and_complete(model,N,data,reps,mu,seconds,nodes,solver_nodes):
    assert len(reps)==len(set(reps))
    payload=f'{model.m} {len(data["labels"])} {model.d}\n'
    payload+=' '.join(map(str,data['distances']))+'\n'
    payload+=' '.join(map(str,data['Gamma']))+'\n'
    payload+=str(len(reps))+'\n'+' '.join(map(str,reps))+'\n'
    cmd=[str(BINARY),'--gamma','--avoid-radical','--first-orbits']
    if seconds is not None:cmd+=['--seconds',str(seconds)]
    if nodes is not None:cmd+=['--max-nodes',str(nodes)]
    tally=Counter();evidence=[];found=None
    with tempfile.TemporaryFile('w+') as err:
        proc=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True)
        try:
            proc.stdin.write(payload);proc.stdin.close()
            for line in proc.stdout:
                fields=line.split();assert fields[0]=='LEAF' and len(fields)==model.m
                xs=[lift_coordinate(int(v),data['units']) for v in fields[1:]]
                assert z_condition(xs,model.Gamma,model.b)
                answer=solve_extension(model,xs,node_limit=solver_nodes,minimum_stabilizer_weight=mu)
                assert answer['status'] in ('FOUND','EXCLUDED','LIMIT')
                tally[answer['status']+':'+answer.get('reason','')]+=1
                evidence.append(dict(X=xs,result=answer))
                if answer['status']=='FOUND':
                    assert model.criterion(answer['model_columns'])
                    S,D=model.construct(answer['model_columns'])
                    assert all(symp(a,b,model.n)==0 for a in S for b in S)
                    assert min(weight(w,model.n) for w in span(S) if w)>=mu
                    assert min(weight(w,model.n) for w in span(D) if w)>=model.d
                    found=dict(S=S,D=D,answer=answer);break
            if found:
                proc.terminate();proc.wait();stats=dict(status='FOUND_EARLY',leaves=len(evidence))
            else:
                assert proc.wait()==0
                err.seek(0);stats=json.loads(err.read())
                assert stats['status'] in ('COMPLETE','LIMIT') and stats['leaves']==len(evidence)
        finally:
            if proc.poll() is None:proc.terminate();proc.wait()
    exhausted=stats['status']=='COMPLETE' and all(e['result']['status']=='EXCLUDED' for e in evidence)
    return dict(status='FOUND' if found else ('EXCLUDED_FIXED_N' if exhausted else 'LIMIT'),
                N=list(N),z=stats,tally=dict(tally),leaves=evidence,found=found,
                minimum_stabilizer_weight=mu,first_representatives=list(reps))


def controls():
    tested=[]
    for n,k,d,m,C,mu in [(5,0,3,3,[],3),(7,1,3,4,[7],3)]:
        model=ExtensionModel(n,k,d,m,[((1<<m)-1)<<m],C)
        data=generic_quotient(model,[])
        reps=[v for v,w in enumerate(data['distances']) if w>=d-1]
        result=enumerate_and_complete(model,[],data,reps,mu,None,None,None)
        assert result['status']=='FOUND'
        tested.append(dict(parameters=[n,k,d],S=result['found']['S'],minimum_stabilizer_weight=mu))
    return tested


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--controls',action='store_true')
    parser.add_argument('--max-tasks',type=int,default=100)
    parser.add_argument('--seconds',type=float,default=1)
    parser.add_argument('--nodes',type=int,default=250000)
    parser.add_argument('--solver-nodes',type=int,default=2000)
    args=parser.parse_args()
    assert __debug__
    expected=fingerprints()
    if args.controls:
        report=dict(status='PASS',positive_controls=controls(),hashes=expected)
        atomic_json(HERE/'verification.json',report)
        print(json.dumps({k:v for k,v in report.items() if k!='hashes'}));return
    verified=json.loads((HERE/'verification.json').read_text())
    assert verified['status']=='PASS' and verified['hashes']==expected
    plan=json.loads((HERE/'plane_plan.json').read_text())
    assert plan['status']=='PASS' and plan['minimum_stabilizer_weight']==6
    assert plan['source_sha256']==sha(SOURCE) and plan['preparation_sha256']==sha(HERE/'prepare.py')
    (HERE/'runs').mkdir(exist_ok=True)
    count=0
    for row in plan['records']:
        model=target(tuple(row['C']))
        assert model.hull==5 and model.C_distance>=6
        for orbit_index,orbit in enumerate(row['orbits']):
            path=HERE/'runs'/f'{row["cls"]}_{orbit_index}.json'
            if path.exists():
                previous=json.loads(path.read_text())
                assert previous['hashes']==expected
                if previous['status']=='EXCLUDED_FIXED_N':continue
            start=time.monotonic();N=orbit['representative']
            data=quotient_input(model,N)
            orbit_data=syndrome_orbits(model,data,N)
            result=enumerate_and_complete(model,N,data,orbit_data['first_representatives'],6,args.seconds,args.nodes,args.solver_nodes)
            result.update(cls=row['cls'],orbit_index=orbit_index,hashes=expected,
                          orbit_data=orbit_data,seconds=time.monotonic()-start,
                          scope='One fixed N orbit in the classified recentered hull-five cell; not a target exclusion')
            atomic_json(path,result)
            print(json.dumps({k:result[k] for k in ('cls','orbit_index','status','seconds','z','tally')}),flush=True)
            count+=1
            if result['status']=='FOUND' or count>=args.max_tasks:return


if __name__=='__main__':main()
