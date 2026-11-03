"""Exact cap oracle, complete leaf comparisons and physical positive controls."""
import itertools
import json
import random
import subprocess
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel, weight,span
from fixed_suffix_search import quotient_input
from signed_extension import solve_extension
from mixed_frames.stage3.common import atomic_json,sha
from rankone_light_budget.common import BINARY,OLD,cap,bound,physical_data,payload


def call(model,data,reps,light=None,budget=None,*,binary=BINARY,xdom=False,limit=None):
    cmd=[str(binary),'--gamma','--avoid-radical','--first-orbits']
    if xdom:cmd+=['--xdom']
    if light is not None:cmd+=['--light-budget',str(budget),'--light-cut',str(light['cut'])]
    if limit is not None:cmd+=['--max-nodes',str(limit)]
    p=subprocess.run(cmd,input=payload(model,data,reps,light),text=True,capture_output=True,check=True,timeout=30)
    stats=json.loads(p.stderr);xs=[tuple(map(int,l.split()[1:])) for l in p.stdout.splitlines()]
    assert stats['leaves']==len(xs)==len(set(xs))
    return stats,xs


def cap_oracle():
    checked=0
    for m in range(3,9):
        for j in range(m):
            values=[-1]*(1<<j)
            for x in range(1<<m):
                if x.bit_count()%2:continue
                for z in range(1<<m):
                    w=x|(z<<m)
                    minimum=min(weight(w,m),weight(w^(((1<<m)-1)<<m),m))
                    low=x&((1<<j)-1)
                    values[low]=max(values[low],minimum)
            assert values==[cap(m,j,low.bit_count()) for low in range(1<<j)]
            checked+=len(values)
    return checked


def main():
    report=dict(status='PASS',cap_patterns=cap_oracle(),complete_comparisons=[],controls=[])
    # Nonempty complete trees with two different physical Gram radicals.
    fixtures=[ExtensionModel(5,0,3,3,[7<<3],[]),ExtensionModel(7,1,3,4,[15<<4],[7])]
    for model in fixtures:
        data=quotient_input(model,[])
        reps=[v for v,d in enumerate(data['distances']) if d>=model.d-1]
        for xdom in (False,True):
            old,L=call(model,data,reps,binary=OLD,xdom=xdom)
            disabled,D=call(model,data,reps,xdom=xdom)
            assert old['status']==disabled['status']=='COMPLETE' and L==D and old['nodes']==disabled['nodes']
            assert L
            for cut in (3,4,5):
                light=physical_data(model,[],data,cut)
                for budget in (0,1,2,10,100):
                    new,K=call(model,data,reps,light,budget,xdom=xdom)
                    assert new['status']=='COMPLETE'
                    # At every earlier depth, the cap can only decrease.
                    # Thus final-bound filtering must equal recursive pruning.
                    expected=[xs for xs in L if bound(xs,light['entries'],model.m,cut)<=budget]
                    assert K==expected,(model.n,cut,budget,xdom,K,expected)
                    report['complete_comparisons'].append(dict(n=model.n,cut=cut,budget=budget,xdom=xdom,plain=len(L),retained=len(K),prunes=new.get('light_pruned',0)))
        # Every physical positive Y completion is retained at its actual
        # measured global light budget, in both native domain-filter lanes.
        positives=0;active_rejections=0
        for xs in L:
            answer=solve_extension(model,list(xs),minimum_stabilizer_weight=model.d)
            if answer['status']!='FOUND':continue
            S,D=model.construct(answer['model_columns'])
            assert min(weight(w,model.n) for w in span(D) if w)>=model.d
            for cut in (model.d+1,model.m+2):
                actual=sum(weight(w,model.n)<cut for w in span(S) if w)
                light=physical_data(model,[],data,cut)
                budget=actual-light['baseline'];assert budget>=0
                for xdom in (False,True):
                    st,K=call(model,data,reps,light,budget,xdom=xdom)
                    assert st['status']=='COMPLETE' and xs in K
                    st,K=call(model,data,reps,light,0,xdom=xdom)
                    active_rejections+=xs not in K
            positives+=1
        assert positives and active_rejections
        report['controls'].append(dict(parameters=[model.n,model.k,model.d],positive_maps=positives,active_underbudget_rejections=active_rejections))
    # No node-limited traversal can be promoted to exhaustion.
    model=fixtures[0];data=quotient_input(model,[]);reps=[v for v,d in enumerate(data['distances']) if d>=2]
    light=physical_data(model,[],data,4)
    for limit in (0,1):
        s,_=call(model,data,reps,light,100,limit=limit)
        assert s['status']=='LIMIT'
    report['hashes']={str(p.relative_to(ROOT)):sha(p) for p in [BINARY,HERE/'z_enum_light.c',HERE/'common.py',Path(__file__),OLD,ROOT/'model.py',ROOT/'signed_extension.py']}
    atomic_json(HERE/'verification.json',report)
    print(json.dumps(dict(status='PASS',cap_patterns=report['cap_patterns'],complete_comparisons=len(report['complete_comparisons']),controls=report['controls'])),flush=True)


if __name__=='__main__':main()
