"""Physically verified local-Clifford automorphisms and syndrome orbits."""
from pynauty import Graph,autgrp
from model import span,symp,lift_coordinate,coordinates


def automorphisms(C,n,N=()):
    words=span(C);W=len(words);Ns=set(span(N));count=W+4*n
    adj={i:set() for i in range(count)}
    def edge(a,b):adj[a].add(b);adj[b].add(a)
    for row,w in enumerate(words):
        for j in range(n):
            letter=((w>>j)&1)|(((w>>(n+j))&1)<<1)
            if letter:edge(row,W+3*j+letter-1)
    for j in range(n):
        for letter in range(3):edge(W+3*n+j,W+3*j+letter)
    colors=[{i for i,w in enumerate(words) if w in Ns},
            {i for i,w in enumerate(words) if w not in Ns},
            set(range(W,W+3*n)),set(range(W+3*n,count))]
    colors=[c for c in colors if c]
    result=autgrp(Graph(count,adjacency_dict=adj,vertex_coloring=colors))
    maps=[]
    for permutation in result[0]:
        sites=[];letters=[]
        for j in range(n):
            targets=[permutation[W+3*j+s]-W for s in range(3)]
            assert len({v//3 for v in targets})==1
            sites.append(targets[0]//3);letters.append([0]+[v%3+1 for v in targets])
            assert sorted(letters[-1])==[0,1,2,3]
        assert sorted(sites)==list(range(n))
        def apply(w,sites=sites,letters=letters):
            out=0
            for j in range(n):
                v=letters[j][((w>>j)&1)|(((w>>(n+j))&1)<<1)]
                out|=(v&1)<<sites[j];out|=(v>>1)<<(n+sites[j])
            return out
        assert {apply(w) for w in words}==set(words)
        assert {apply(w) for w in Ns}==Ns
        units=[1<<i for i in range(2*n)]
        assert all(symp(apply(a),apply(b),n)==symp(a,b,n) for a in units for b in units)
        maps.append(apply)
    return maps,round(result[1]*10**result[2])


def syndrome_orbits(model,data,N):
    maps,order=automorphisms(model.C,model.f,N)
    labels=data['labels'];units=data['units'];size=len(labels)
    parent=list(range(size));transformations=[]
    def find(v):
        while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
        return v
    for apply in maps:
        image=[]
        for old in labels:
            new=model.syndrome(apply(model.representatives[old]))
            image.append(coordinates(new,units))
        assert len(set(image))==size
        columns=[image[1<<i] for i in range(len(units))]
        assert all(image[v]==lift_coordinate(v,columns) for v in range(size))
        assert all(data['distances'][v]==data['distances'][image[v]] for v in range(size))
        # T Gamma T^T = Gamma, with T acting on syndrome functionals.
        rows=[sum(((columns[j]>>i)&1)<<j for j in range(len(units))) for i in range(len(units))]
        TG=[lift_coordinate(row,data['Gamma']) for row in rows]
        transformed=[sum(((TG[i]&rows[j]).bit_count()%2)<<j for j in range(len(units)))
                     for i in range(len(units))]
        assert transformed==data['Gamma']
        for v,w in enumerate(image):
            a,b=find(v),find(w)
            if a!=b:parent[max(a,b)]=min(a,b)
        transformations.append(columns)
    high=[i for i,d in enumerate(data['distances']) if d>=5]
    reps=sorted({find(v) for v in high})
    assert all(data['distances'][v]>=5 for v in reps)
    return dict(first_representatives=reps,high_labels=len(high),first_orbits=len(reps),
                automorphism_order=order,verified_syndrome_generators=transformations,
                scope='C,N,Pauli weights, symplectic form, coset distance and syndrome linearity preserved')
