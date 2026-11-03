"""Physical suffix automorphisms, dual syndrome actions, and affine normal forms."""
from collections import deque
from model import coordinates, lift_coordinate, reduced_basis, span, symp, weight
from pynauty import Graph, autgrp


def physical_apply(w, sites, letters, n):
    out=0
    for j in range(n):
        a=((w>>j)&1)|(((w>>(n+j))&1)<<1)
        s=letters[j][a]
        out|=(s&1)<<sites[j];out|=(s>>1)<<(n+sites[j])
    return out


def automorphisms(m):
    n=m.f;words=span(m.C);W=len(words)
    adj={i:set() for i in range(W+3*n)}
    def edge(a,b):adj[a].add(b);adj[b].add(a)
    for row,w in enumerate(words):
        for j in range(n):
            s=((w>>j)&1)|(((w>>(n+j))&1)<<1)
            if s:edge(row,W+3*j+s-1)
    for j in range(n):
        for a in range(3):
            for b in range(a):edge(W+3*j+a,W+3*j+b)
    autos=autgrp(Graph(W+3*n,adjacency_dict=adj,vertex_coloring=[set(range(W)),set(range(W,W+3*n))]))
    gens=[]
    for perm in autos[0]:
        sites=[];letters=[]
        for j in range(n):
            targets=[perm[W+3*j+s]-W for s in range(3)]
            assert len({v//3 for v in targets})==1
            sites.append(targets[0]//3);letters.append([0]+[v%3+1 for v in targets])
        assert sorted(sites)==list(range(n)) and all(sorted(x)==[0,1,2,3] for x in letters)
        apply=lambda w:physical_apply(w,sites,letters,n)
        assert {apply(w) for w in words}==set(words)
        T=[coordinates(apply(w),m.P) for w in m.P]
        assert len(reduced_basis(T))==m.b
        inverse=[coordinates(1<<i,T) for i in range(m.b)]
        S=[sum(((inverse[i]>>j)&1)<<i for i in range(m.b)) for j in range(m.b)]
        assert all(((T[i]&S[j]).bit_count()&1)==int(i==j) for i in range(m.b) for j in range(m.b))
        assert all(symp(apply(a),apply(b),n)==symp(a,b,n) for a in m.P for b in m.P)
        wt=[weight(lift_coordinate(a,m.P),n) for a in range(1<<m.b)] # span uses binary coefficient order
        for a in range(1<<m.b):
            ta=lift_coordinate(a,T);ss=lift_coordinate(a,S)
            assert wt[a]==wt[ta] and m.distances[a]==m.distances[ss]
            assert lift_coordinate(lift_coordinate(a,m.Gamma),S)==lift_coordinate(ta,m.Gamma)
            # Covers every syndrome, even outside im(Gamma) for singular forms.
            assert m.syndrome(apply(m.representatives[a]))==ss
        gens.append(dict(sites=sites,letters=letters,coefficient_columns=T,syndrome_columns=S))
    return dict(generators=gens,graph_group_order=round(autos[1]*10**autos[2]),hull=m.hull)


def coverage(eligible, maps):
    """Forward BFS certificates: parent --generator--> member from each minimum."""
    remaining=set(eligible);orbits=[]
    while remaining:
        root=min(remaining);parent={root:None};todo=deque([root])
        while todo:
            a=todo.popleft()
            for k,T in enumerate(maps):
                b=lift_coordinate(a,T)
                assert b in remaining
                if b not in parent:parent[b]=[a,k];todo.append(b)
        remaining-=parent.keys()
        orbits.append(dict(representative=root,members=sorted(parent),witnesses={str(v):parent[v] for v in sorted(parent)}))
    return orbits


def affine_points(xs):
    return [0,*xs,lift_coordinate((1<<len(xs))-1,xs)]


def anchor_normal(xs):
    pts=affine_points(xs)
    best=min(sorted(v^anchor for v in pts) for anchor in pts)
    assert len(set(best))==len(xs)+2
    assert best[-1]==lift_coordinate((1<<len(xs))-1,best[1:-1])
    return best[1:-1]


def joint_normal(xs, cert):
    """Construct a safe representative without enumerating the full suffix group.

    Move a smallest edge-orbit minimum to its representative, then normalize
    anchors. This realizes the first component of the joint lex minimum;
    the remaining suffix stabilizer ties need not be removed for coverage.
    """
    membership={v:o for o in cert['orbits'] for v in o['members']}
    pts=affine_points(xs)
    edges={a^b for i,a in enumerate(pts) for b in pts[:i]}
    chosen=min(edges,key=lambda v:(membership[v]['representative'],v))
    orbit=membership[chosen];path=[];v=chosen
    while orbit['witnesses'][str(v)] is not None:
        parent,k=orbit['witnesses'][str(v)];path.append(k);v=parent
    out=list(xs)
    for k in path:
        T=cert['inverse_syndrome_columns'][k]
        out=[lift_coordinate(a,T) for a in out]
    out=anchor_normal(out)
    assert out[0]==orbit['representative']
    return out


def validate_certificate(cert,aut):
    maps=[g['syndrome_columns'] for g in aut['generators']]
    seen=set()
    for o in cert['orbits']:
        root=o['representative'];members=set(o['members'])
        assert root==min(members) and not seen.intersection(members)
        assert set(o['witnesses'])=={str(v) for v in members}
        children={v:[] for v in members}
        for v in members:
            edge=o['witnesses'][str(v)]
            if v==root:assert edge is None
            else:
                parent,k=edge;assert parent in members and 0<=k<len(maps)
                assert lift_coordinate(parent,maps[k])==v
                children[parent].append(v)
            for S in maps:assert lift_coordinate(v,S) in members
        reached={root};todo=[root]
        while todo:
            for v in children[todo.pop()]:
                assert v not in reached;reached.add(v);todo.append(v)
        assert reached==members
        seen|=members
    assert seen==set(cert['eligible']) and len(seen)==len(cert['eligible'])
    assert cert['first_representatives']==[o['representative'] for o in cert['orbits']]
    assert cert['inverse_syndrome_columns']==[[coordinates(1<<j,S) for j in range(len(S))] for S in maps]
