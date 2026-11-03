"""Independent cross-check of suffix9_dim5.json: Aut-orbit-reduced, labelled-mass-checked
lengthening census at distance 6 (coset min >=5) directly from the empty length-5 code to
(9,2^5,>=6). Does not use cell_8_5.json or the 187-form lift. Checkpoint per cell."""
import json,time
from fractions import Fraction
from math import factorial
from pathlib import Path
import census8_5 as m
OUT=Path(__file__).resolve().parent
def deep(g,n):
    label,lift,q=m.quotient_labels(g,n);cov={label(e) for e in m.light_errors(n,4)}
    return label,lift,[a for a in range(1<<q) if a not in cov]
m.deep_cosets=deep
def ext_any(g,r,n):
    return (1,[list(g)]) if r==0 else m.extend(g,r,n)
cells={(5,0):[([],m.aut_order(m.graph([],5)))]};log=[];maxn,tk=9,5
for n in range(6,maxn+1):
    for k in range(max(0,tk-2*(maxn-n)),tk+1):
        t0=time.time();classes={};lhs=Fraction(0)
        for r in (0,1,2):
            for g0,aut in cells.get((n-1,k-r),[]):
                if not g0 and r>0:  # zero parent: Aut generators fine, handled by extend
                    pass
                e,gs=ext_any(g0,r,n-1);lhs+=Fraction(factorial(n-1)*6**(n-1),aut)*e
                for g in gs:
                    G=m.graph(g,n);c=m.certificate(G)
                    if c not in classes:classes[c]=(g,m.aut_order(G))
        found=list(classes.values())
        for g,a in found:
            s=m.span(g);assert len(set(s))==1<<k and (k==0 or min(m.wt(v,n) for v in s if v)>=6)
        rhs=sum(Fraction(factorial(n)*6**n,a) for _,a in found)
        e=dict(cell=[n,k],classes=len(found),mass_equal=lhs==rhs,mass=str(rhs),seconds=round(time.time()-t0,1))
        print(json.dumps(e),flush=True);log.append(e);assert lhs==rhs;cells[(n,k)]=found
        (OUT/'census_d6_orbit.json').write_text(json.dumps(dict(status='RUNNING',cells=log))+'\n')
mine={m.certificate(m.graph(g,9)) for g,_ in cells[(9,5)]}
s=json.loads((OUT/'suffix9_dim5.json').read_text())
inter=lambda v,n:sum((((v>>i)&1)<<(2*i))|(((v>>(n+i))&1)<<(2*i+1)) for i in range(n))
theirs={m.certificate(m.graph([inter(v,9) for v in c['representative']['generators']],9)) for c in s['classes']}
res=dict(status='COMPLETE',cells=log,classes_9_5=len(mine),lift_classes=len(theirs),certificate_sets_equal=mine==theirs)
(OUT/'census_d6_orbit.json').write_text(json.dumps(res,indent=1)+'\n')
print(json.dumps({k:res[k] for k in ['classes_9_5','lift_classes','certificate_sets_equal']}))
