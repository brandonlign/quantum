"""Complete small augmented-prefix cells for the m8,r1 gauge reduction."""
from pathlib import Path
import json,hashlib,sys
from collections import Counter
from pynauty import Graph,certificate
from classify_suffix9 import maps,lift
from model import span,weight,gram,reduced_basis,quotient_basis,lift_coordinate
BASE=Path(__file__).resolve().parent
DATA=BASE.parents[1]/'release16/census/cells.json'


def cert(g,n):
    ws=span(g);W=len(ws);adj={i:set() for i in range(W+3*n)}
    def edge(a,b):adj[a].add(b);adj[b].add(a)
    for row,w in enumerate(ws):
        for j in range(n):
            s=((w>>j)&1)|(((w>>(n+j))&1)<<1)
            if s:edge(row,W+3*j+s-1)
    for j in range(n):
        for a in range(3):
            for b in range(a+1,3):edge(W+3*j+a,W+3*j+b)
    return certificate(Graph(W+3*n,adjacency_dict=adj,vertex_coloring=[set(range(W)),set(range(W,W+3*n))])).hex()


def from_inter(w,n):
    x=sum(((w>>(2*i))&1)<<i for i in range(n));z=sum(((w>>(2*i+1))&1)<<i for i in range(n))
    return x|(z<<n)


def main():
    cells=json.loads(DATA.read_text());groups={}
    for dim,parents in [(3,json.loads((BASE/'opus_parallel/suffix9_dimension5/cells_d5_7_3.json').read_text())['7,3']['xz']),(4,cells['7,4'])]:
        records={};reps=0
        for pi,g in enumerate(parents):
            ws=[0]
            for v in g:ws += [w^v for w in ws]
            assert len(set(ws))==1<<dim and min(weight(w,7) for w in ws if w)>=5
            bad=[i for i,w in enumerate(ws) if weight(w,7)==5]
            for a,b in maps(dim):
                if not all((a&w).bit_count()%2 or (b&w).bit_count()%2 for w in bad):continue
                R=lift(g,a,b,7);code=span(R)
                assert min(weight(w,8) for w in code if w)>=6
                if len(reduced_basis(gram(R,8))) or not any(weight(w,8)==8 for w in code):continue
                reps+=1;key=cert(R,8)
                if key not in records:
                    E=quotient_basis(R,8);profile=Counter(min(weight(lift_coordinate(e,E)^r,8) for r in code) for e in range(1<<len(E)))
                    records[key]=dict(generators=R,weight_distribution=dict(Counter(weight(w,8) for w in code)),
                        quotient_dimension=len(E),cost_profile=dict(profile),max_cost=max(profile))
        groups[str(dim)]=dict(status='COMPLETE',parent_classes=len(parents),canonical_extension_representatives=reps,
                              inequivalent_classes=len(records),classes=list(records.values()))
    p=255<<8
    R2=[[p,15|(48<<8)],[p,63],[p,63|(64<<8)],[p,255]]
    records=[]
    for R in R2:
        code=span(R);E=quotient_basis(R,8);profile=Counter(min(weight(lift_coordinate(e,E)^r,8) for r in code) for e in range(1<<len(E)))
        records.append(dict(generators=R,weight_distribution=dict(Counter(weight(w,8) for w in code)),
                            quotient_dimension=len(E),cost_profile=dict(profile),max_cost=max(profile)))
    groups['2']=dict(status='COMPLETE_DISTINGUISHED_R_FORMS',inequivalent_classes=4,classes=records)
    result=dict(status='COMPLETE_PREFIX_CANDIDATES',groups=groups,
        scope='Prefix local-code candidates for a proposed maximal-isotropic kernel gauge reduction; '
              'does not by itself establish an original code extension or parameter exclusion.',
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),DATA,
            BASE/'opus_parallel/suffix9_dimension5/cells_d5_7_3.json',BASE/'opus_parallel/suffix9_dimension5/verify73_types.json']})
    (BASE/'gauge_prefix_candidates.json').write_text(json.dumps(result,indent=2)+'\n')
    for k,v in groups.items():print(k,{a:b for a,b in v.items() if a!='classes'},[r['cost_profile'] for r in v['classes']])
if __name__=='__main__':main()
