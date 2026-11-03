"""Kernel-zero state-base tables for the recentered hull-five rank-one cell."""
from functools import lru_cache
from pathlib import Path
from model import (ExtensionModel,coordinates,gram,lift_coordinate,reduced_basis,
                   span,weight,symp)
from mixed_frames.frames import Frame,adapted_frame,xor_table
from mixed_frames.stage3.common import atomic_json,sha

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
BINARY=ROOT/'mixed_frames/stage3/enumerate_z'
SOURCE=ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
PREFIX=ROOT/'gauge_prefix_candidates.json'


def normalize_word(w,p,m):
    out=0
    for i in range(m):
        letter=((w>>i)&1)|(((w>>(m+i))&1)<<1)
        distinguished=((p>>i)&1)|(((p>>(m+i))&1)<<1)
        assert distinguished
        table={1:[0,2,1,3],2:[0,1,2,3],3:[0,1,3,2]}[distinguished]
        v=table[letter]
        out|=(v&1)<<i;out|=(v>>1)<<(m+i)
    return out


def cheap_frame(m,R,name):
    initial=adapted_frame(m,R,name)
    p=initial.p;rs=span(initial.R)
    choices=sorted(range(1,1<<p),key=lambda z:(initial.costs[z],lift_coordinate(z,initial.E[:p])))
    columns=[]
    for z in choices:
        if len(reduced_basis(columns+[z]))>len(columns):columns.append(z)
        if len(columns)==p:break
    inverse=[coordinates(1<<i,columns) for i in range(p)]
    paired=[sum(((inverse[i]>>j)&1)<<i for i in range(p)) for j in range(p)]
    E=[lift_coordinate(z,initial.E[:p]) for z in columns]+[lift_coordinate(x,initial.E[p:]) for x in paired]
    assert gram(E,m)==[1<<(p+i) for i in range(p)]+[1<<i for i in range(p)]
    costs=[min(weight(w^r,m) for r in rs) for w in xor_table(E)]
    return Frame(m,initial.R,E,costs,name)


def prefix_candidates(records):
    candidates=[]
    for index,row in enumerate(records):
        H=reduced_basis(row['generators'])
        assert len(H)==4 and min(weight(w,8) for w in span(H) if w)>=6
        for distinguished in sorted(w for w in span(H) if weight(w,8)==8):
            R=reduced_basis([normalize_word(w,distinguished,8) for w in H])
            assert normalize_word(distinguished,distinguished,8)==255<<8
            assert min(weight(w,8) for w in span(R) if w)>=6
            assert all(symp(a,b,8)==0 for a in R for b in R)
            assert [w for w in span(R) if w and not w&255]==[255<<8]
            frame=cheap_frame(8,R,f'h5_gauge_{len(candidates)}')
            assert frame.p==4 and max(frame.cost[0])<6
            candidates.append(dict(prefix_index=len(candidates),class_index=index,
                                   distinguished_word=distinguished,R=R,E=frame.E))
    assert len(candidates)==16
    return candidates


def quotient(model):
    # The state-base shortened C is exactly the radical of P=C-perp.
    C=model.C;h=len(C)
    assert all(symp(c,p,model.f)==0 for c in C for p in model.P)
    full=[coordinates(c,model.P) for c in C];complement=[]
    for i in range(model.b):
        v=1<<i
        if len(reduced_basis(full+[v]))>len(full):full.append(v);complement.append(v)
    a=len(complement)
    assert a==model.a and len(full)==model.b
    Q=[lift_coordinate(v,model.P) for v in complement]
    G=gram(Q,model.f);assert len(reduced_basis(G))==a
    units=[sum(((coordinates(1<<i,full)>>(h+j))&1)<<i for i in range(model.b)) for j in range(a)]
    labels=xor_table(units)
    weights=[min(weight(w^c,model.f) for c in span(C)) for w in xor_table(Q)]
    distances=[model.distances[s] for s in labels]
    assert len(set(labels))==1<<a
    # Independently verify the full coefficient Gram factors through Q.
    for i in range(model.b):
        pi=sum(((u>>i)&1)<<j for j,u in enumerate(units))
        for j in range(model.b):
            pj=sum(((u>>j)&1)<<k for k,u in enumerate(units))
            assert ((model.Gamma[i]>>j)&1)==((pi&lift_coordinate(pj,G)).bit_count()%2)
    for c,w in enumerate(weights):
        assert w==distances[lift_coordinate(c,G)]
    return dict(Gamma=G,units=units,labels=labels,distances=distances,weights=weights,Q=Q)


@lru_cache(8)
def target(C,R,name):
    frame=cheap_frame(8,list(R),name)
    model=ExtensionModel(17,0,6,8,frame.R,list(C))
    assert model.a==8 and model.b==13 and model.hull==5
    return model,frame


def suffix_generators(model,data):
    from m8r1.suffix_symmetries import automorphisms
    maps,order=automorphisms(model.C,model.f,model.C)
    out=[];a=len(data['Gamma']);size=1<<a
    for apply in maps:
        image=[coordinates(model.syndrome(apply(model.representatives[s])),data['units']) for s in data['labels']]
        columns=[image[1<<i] for i in range(a)]
        assert len(set(image))==size
        assert all(image[v]==lift_coordinate(v,columns) for v in range(size))
        assert all(data['distances'][v]==data['distances'][image[v]] for v in range(size))
        rows=[sum(((columns[j]>>i)&1)<<j for j in range(a)) for i in range(a)]
        transformed=[sum(((lift_coordinate(rows[i],data['Gamma'])&rows[j]).bit_count()%2)<<j for j in range(a)) for i in range(a)]
        assert transformed==data['Gamma']
        out.append(columns)
    return out,order


def first_orbits(eligible,generators):
    remaining=set(eligible);out=[]
    while remaining:
        root=min(remaining);seen={root};todo=[root]
        while todo:
            v=todo.pop()
            for columns in generators:
                w=lift_coordinate(v,columns);assert w in remaining
                if w not in seen:seen.add(w);todo.append(w)
        remaining-=seen
        out.append(dict(representative=root,members=sorted(seen)))
    return out


def input_text(model,frame,data,mu=None):
    if mu is None:mu=model.d
    arrays=[[frame.p,len(data['Gamma']),model.d,mu,1,1],frame.costs,
            data['distances'],data['Gamma'],data['weights']]
    return '\n'.join(' '.join(map(str,row)) for row in arrays)+'\n'


def hashes():
    paths=[Path(__file__),HERE/'prepare.py',HERE/'verify.py',HERE/'run.py',HERE/'PROOF.md',
           SOURCE,PREFIX,BINARY,ROOT/'mixed_frames/stage3/enumerate_z.cpp',
           ROOT/'mixed_frames/frames.py',ROOT/'mixed_frames/fixed_z.py',ROOT/'model.py',
           ROOT/'r1_column/Q4_KERNEL.md',ROOT/'CRITICAL_PREFIX_SPECTRUM.md',ROOT/'gauge_fix.py']
    return {str(p.relative_to(ROOT)):sha(p) for p in paths}
