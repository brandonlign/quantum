#!/usr/bin/env python3
"""Fresh finite orbit partitions for the full-plane rank-two case.

Uses pynauty's full physical graph automorphism generators. Enumerates every
subspace explicitly, checks the partition mass, and separately computes all
168 radical identifications. No production search modules are imported.
"""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import json,time
import pynauty
from check_local_classification import span,weight,symp
HERE=Path(__file__).resolve().parent

def automorphisms(rows,n):
    words=span(rows);W=len(words);end=W+3*n;adj={i:set() for i in range(end+n)}
    def edge(a,b):adj[a].add(b);adj[b].add(a)
    for j,w in enumerate(words):
        for i in range(n):
            s=((w>>i)&1)|(((w>>(n+i))&1)<<1)
            if s:edge(j,W+3*i+s-1)
    for i in range(n):
        for s in range(3):edge(end+i,W+3*i+s)
    g=pynauty.Graph(end+n,adjacency_dict=adj,vertex_coloring=[set(range(W)),set(range(W,end)),set(range(end,end+n))])
    perms,_,_,_,_=pynauty.autgrp(g)
    maps=[]
    for perm in perms:
        action={w:words[perm[j]] for j,w in enumerate(words)}
        assert action[0]==0
        assert all(action[u^v]==action[u]^action[v] for u in words for v in rows)
        assert all(weight(w,n)==weight(action[w],n) for w in words)
        assert all(symp(u,v,n)==symp(action[u],action[v],n) for u in rows for v in rows)
        maps.append(action)
    return maps

def partition(objects,actions):
    lookup={x:i for i,x in enumerate(objects)};assert len(lookup)==len(objects)
    parent=list(range(len(objects)))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for action in actions:
        images=[lookup[action(x)] for x in objects];assert sorted(images)==list(range(len(objects)))
        for i,j in enumerate(images):
            a,b=root(i),root(j)
            if a!=b:parent[b]=a
    out={}
    for i,obj in enumerate(objects):out.setdefault(root(i),[]).append(obj)
    parts=list(out.values());assert sum(map(len,parts))==len(objects)
    return parts

def main():
    (Path(__file__).resolve().parents[1]/'audit').mkdir(exist_ok=True)
    assert __debug__
    start=time.monotonic();rows=json.loads((HERE/'data/ranktwo-shortened-codes.json').read_text())
    counts=Counter();records=[]
    for row in rows:
        C=row['C'];h=row['hull'];q=row['q'];words=span(C)
        assert len(C)==6 and min(weight(w,9) for w in words[1:])>=6
        radical=[w for w in words if all(symp(w,v,9)==0 for v in C)]
        assert len(radical)==1<<h
        if not q:continue
        nonzero=radical[1:]
        if q==1:spaces=[(w,) for w in nonzero];assert len(spaces)==3
        else:
            spaces=sorted({tuple(span([a,b,c])[1:]) for a,b,c in combinations(nonzero,3) if a^b!=c})
            assert len(spaces)==1395
        maps=automorphisms(C,9)
        parts=partition(spaces,[lambda N,m=m:tuple(sorted(m[w] for w in N)) for m in maps])
        tally=Counter()
        for orbit in parts:
            flags={sum(weight(w,9) in (6,7) for w in N)>2 for N in orbit}
            assert len(flags)==1
            tally['screened' if flags.pop() else 'remaining']+=1
        for k,v in tally.items():counts[f'q{q}_{k}']+=v
        records.append({'class':row['cls'],'q':q,'all_subspaces':len(spaces),'orbits':len(parts),'orbit_sizes':sorted(map(len,parts)),'counts':dict(tally)})
    assert counts['q1_remaining']==63 and counts['q3_remaining']==164,counts
    local=json.loads((HERE.parent/'audit/review-local-classification.json').read_text())
    H=local['full_plane_dimension5_representatives'][0];full=[w for w in span(H) if weight(w,8)==8];R={0,*full}
    quotient={};representatives=[]
    for w in span(H):
        if w in quotient:continue
        i=len(representatives);representatives.append(w)
        for r in R:quotient[w^r]=i
    assert len(representatives)==8
    # Identify quotient addition independently of its numerical labeling.
    add=lambda a,b:quotient[representatives[a]^representatives[b]]
    ids=[t for t in product(range(1,8),repeat=3) if t[0]!=t[1] and t[2] not in (t[0],t[1],add(t[0],t[1]))]
    assert len(ids)==168
    maps=automorphisms(H,8);actions=[]
    for m in maps:
        assert {m[r] for r in R}==R
        table=tuple(quotient[m[w]] for w in representatives)
        actions.append(table)
    identity=tuple(range(8));group={identity};todo=[identity]
    while todo:
        a=todo.pop()
        for b in actions:
            c=tuple(b[a[i]] for i in range(8))
            if c not in group:group.add(c);todo.append(c)
    assert len(group)==21
    idparts=partition(ids,[lambda t,g=g:tuple(g[i] for i in t) for g in actions])
    assert len(idparts)==8 and all(len(p)==21 for p in idparts)
    report={'status':'PASS','scope':'Complete fresh full-plane subspace and radical-identification orbits, conditional on the 56 shortened-code representatives and local classification.',
        'counts':dict(counts),'classes':len(rows),'nonzero_hull_classes':len(records),'records':records,
        'radical_identifications':168,'physical_quotient_group_order':21,'identification_orbits':8,'identification_orbit_sizes':list(map(len,idparts)),
        'seconds':time.monotonic()-start}
    (HERE.parent/'audit/review-full-plane-orbits.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
if __name__=='__main__':main()
