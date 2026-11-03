"""Construct the ordered state basis and verify exact inverse selections."""
import json
from pathlib import Path
import subprocess
from model import (coordinates,dual,embed,lift_coordinate,nullspace,reduced_basis,
                   span,symp,weight)

HERE=Path(__file__).resolve().parent


def ordered_basis(model,frame,columns):
    p=((1<<model.m)-1)<<model.m
    H=[p]
    for r in frame.R:
        if len(reduced_basis(H+[r]))>len(H):H.append(r)
    assert len(H)==len(frame.R)
    B=[embed(h,0,model.m,model.f) for h in H]
    B += [embed(0,c,model.m,model.f) for c in model.C]
    B += [embed(e,model.representatives[s],model.m,model.f) for e,s in zip(frame.E,columns)]
    assert len(B)==model.n and len(reduced_basis(B))==model.n
    assert all(symp(a,b,model.n)==0 for a in B for b in B)
    return B


def select(B,n,m,h,d,cut,budget,seconds=None,nodes=None):
    text=f'{n} {d} {m} {h} {cut} {budget}\n'+' '.join(map(str,B))+'\n'
    cmd=[str(HERE/'inverse')]
    if seconds is not None:cmd+=['--seconds',str(seconds)]
    if nodes is not None:cmd+=['--nodes',str(nodes)]
    p=subprocess.run(cmd,input=text,text=True,capture_output=True,check=True)
    result=json.loads(p.stdout)
    assert result['status'] in ('FOUND','EXCLUDED_BASE','LIMIT')
    if result['status']=='FOUND':
        rows=result['rows'];assert len(rows)==h-1
        assert all((row&((1<<h)-1))==1<<(i+1) for i,row in enumerate(rows))
        S=reduced_basis([lift_coordinate(c,B) for c in nullspace(rows,n)])
        assert len(S)==n-h+1
        assert all(symp(a,b,n)==0 for a in S for b in S)
        D=dual(S,n)
        assert min(weight(w,n) for w in span(D) if w)>=d
        assert sum(weight(w,n)<cut for w in span(S) if w)<=budget
        from gauge_fix import supported_subspace
        R=[embed(((1<<m)-1)<<m,0,m,n-m)]
        assert supported_subspace(S,n,range(m))==reduced_basis(R)
        result.update(generators=S,normalizer=D)
    return result
