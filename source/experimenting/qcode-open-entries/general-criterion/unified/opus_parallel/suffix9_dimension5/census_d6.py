"""Independent check: mass-checked lengthening census at distance 6 directly to (9,2^5,>=6),
using mass_census.py read-only with D replaced by 6 (coset min >=5). Compares the
certificate set with suffix9_dim5.json classes. Checkpoints per cell."""
import json,time
from fractions import Fraction
from math import factorial
from pathlib import Path
OUT=Path(__file__).resolve().parent
SRC=Path('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/mass_census.py')
src=SRC.read_text().split('if __name__ == "__main__":')[0].replace('D = 5  # target','D = 6  # target')
ns={'__file__':str(SRC)};exec(src,ns);assert ns['D']==6
span,wt,graph,aut_order,extensions,certificate=(ns[k] for k in ['span','wt','graph','aut_order','extensions','certificate'])
cells={(5,0):[([],aut_order(graph([],5)))]};log=[]
maxn,tk=9,5
for n in range(6,maxn+1):
    for k in range(max(0,tk-2*(maxn-n)),tk+1):
        t0=time.time();classes={};lhs=Fraction(0)
        for r in (0,1,2):
            for g0,aut in cells.get((n-1,k-r),[]):
                exts,ext=extensions(g0,n-1);lhs+=Fraction(factorial(n-1)*6**(n-1),aut)*ext[r]
                for g in exts[r]:
                    G=graph(g,n);c=certificate(G)
                    if c not in classes:classes[c]=(g,aut_order(G))
        found=list(classes.values())
        for g,a in found:
            s=span(g);assert len(set(s))==1<<k and (k==0 or min(wt(v,n) for v in s if v)>=6)
        rhs=sum(Fraction(factorial(n)*6**n,a) for _,a in found)
        e=dict(cell=[n,k],classes=len(found),mass_equal=lhs==rhs,mass=str(rhs),seconds=round(time.time()-t0,1))
        print(json.dumps(e),flush=True);log.append(e);assert lhs==rhs;cells[(n,k)]=found
        (OUT/'census_d6.json').write_text(json.dumps(dict(status='RUNNING',cells=log))+'\n')
mine={certificate(graph(g,9)) for g,_ in cells[(9,5)]}
s=json.loads((OUT/'suffix9_dim5.json').read_text())
def inter(v,n):return sum((((v>>i)&1)<<(2*i))|(((v>>(n+i))&1)<<(2*i+1)) for i in range(n))
theirs={certificate(graph([inter(v,9) for v in c['representative']['generators']],9)) for c in s['classes']}
res=dict(status='COMPLETE',cells=log,classes_9_5=len(mine),lift_classes=len(theirs),certificate_sets_equal=mine==theirs)
(OUT/'census_d6.json').write_text(json.dumps(res,indent=1)+'\n');print(json.dumps({k:res[k] for k in ['classes_9_5','lift_classes','certificate_sets_equal']}))
