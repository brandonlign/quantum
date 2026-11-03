"""Control: rerun census8_5.extend (Aut-orbit reduced) for cell (8,2^6) from (7,2^5)[r1],(7,2^4)[r2]
((7,2^6) is empty by census85 mass check / release16), compare to release16 cells.json '8,6' (1579)."""
import json
from fractions import Fraction
from math import factorial
from pathlib import Path
import census8_5 as m
OUT=Path(__file__).resolve().parent
ps={1:json.loads((OUT/'cells_d5_7_5.json').read_text())['7,5']['interleaved_with_aut'],
    2:json.loads((OUT/'cells_d5_7_4.json').read_text())['7,4']['interleaved_with_aut']}
classes={};lhs=Fraction(0)
for r,L in ps.items():
    for g,a in L:
        ext,gs=m.extend(g,r);lhs+=Fraction(factorial(7)*6**7,a)*ext
        for h in gs:
            G=m.graph(h,8);c=m.certificate(G)
            if c not in classes:classes[c]=m.aut_order(G)
rhs=sum(Fraction(factorial(8)*6**8,a) for a in classes.values())
rel=json.loads(Path('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/cells.json').read_text())['8,6']
inter=lambda v,n:sum((((v>>i)&1)<<(2*i))|(((v>>(n+i))&1)<<(2*i+1)) for i in range(n))
theirs={m.certificate(m.graph([inter(v,8) for v in g],8)) for g in rel}
res=dict(classes=len(classes),release16=len(theirs),sets_equal=set(classes)==theirs,mass_equal=lhs==rhs,mass=str(rhs))
(OUT/'control_8_6.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
assert res['sets_equal'] and res['mass_equal']
