"""Verified suffix-code automorphisms and safe first-Z orbit representatives."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
from pynauty import Graph,autgrp
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from model import span,coordinates,lift_coordinate,symp
BASE=Path(__file__).resolve().parent


def main():
    inputs=json.loads((BASE.parent/'rank2_z_inputs.json').read_text())['classes'];reports=[]
    for item in inputs:
        d=json.loads(Path(item['table']).read_text());C=d['C_generators'];P=d['P_generators'];n=9;ws=span(C);W=len(ws)
        adj={i:set() for i in range(W+3*n)}
        for row,w in enumerate(ws):
            for j in range(n):
                s=((w>>j)&1)|(((w>>(n+j))&1)<<1)
                if s:adj[row].add(W+3*j+s-1)
        for j in range(n):
            for s in range(3):adj[W+3*j+s].update(W+3*j+t for t in range(3) if t!=s)
        autos=autgrp(Graph(W+3*n,adjacency_dict=adj,vertex_coloring=[set(range(W)),set(range(W,W+3*n))]))
        N=len(d['flags']);parent=list(range(N));maps=[];high=[v for v in range(1,N) if d['flags'][v]&1]
        def find(a):
            while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
            return a
        for perm in autos[0]:
            site=[];letters=[]
            for j in range(n):
                targets=[perm[W+3*j+s-1]-W for s in [1,2,3]]
                assert len({v//3 for v in targets})==1
                site.append(targets[0]//3);letters.append([0]+[v%3+1 for v in targets])
                assert sorted(letters[-1])==[0,1,2,3]
            assert sorted(site)==list(range(n))
            def apply(w):
                out=0
                for j in range(n):
                    s=((w>>j)&1)|(((w>>(n+j))&1)<<1);v=letters[j][s]
                    out|=(v&1)<<site[j];out|=(v>>1)<<(n+site[j])
                return out
            assert {apply(w) for w in ws}==set(ws)
            images=[coordinates(apply(w),P) for w in P]
            assert all(symp(apply(a),apply(b),n)==symp(a,b,n) for a in P for b in P)
            transform=[lift_coordinate(a,images) for a in range(N)]
            assert len(set(transform))==N
            for a,b in enumerate(transform):
                assert d['suffix_weights'][a]==d['suffix_weights'][b]
                assert d['coset_distances'][a]==d['coset_distances'][b]
                assert d['flags'][a]==d['flags'][b]
            for a in high:
                ra,rb=find(a),find(transform[a])
                if ra!=rb:parent[max(ra,rb)]=min(ra,rb)
            maps.append(images)
        groups={}
        for a in high:groups.setdefault(find(a),[]).append(a)
        reps=sorted(min(v) for v in groups.values())
        r=dict(class_index=d['class_index'],high_labels=len(high),first_representatives=reps,
            first_orbits=len(reps),verified_generator_coefficient_maps=maps,
            graph_automorphism_order=round(autos[1]*10**autos[2]),
            table_sha256=hashlib.sha256(Path(item['table']).read_bytes()).hexdigest())
        reports.append(r);print(json.dumps({k:v for k,v in r.items() if k not in ['first_representatives','verified_generator_coefficient_maps']}),flush=True)
    result=dict(status='PASS',classes=reports,
        scope='Weight, coset distance, symplectic form and C preservation checked. '
              'For any sorted7-image Z circuit, an automorphism can move the smallest '
              'orbit minimum to the smallest image, so restricting the first image '
              'to these representatives preserves existence of a full extension.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (BASE/'aut_orbits.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
