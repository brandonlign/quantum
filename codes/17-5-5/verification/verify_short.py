import sys
if sys.flags.optimize:
    raise RuntimeError("Run without -O: this certificate uses assertions")
"""Independent short distance certificate, using only the printed Pauli rows."""
from pathlib import Path
from itertools import combinations, product
from collections import Counter

strings = (Path(__file__).parent / 'generators.txt').read_text().splitlines()
assert len(strings) == 12
labels = {'I': (0,0), 'X': (1,0), 'Z': (0,1), 'Y': (1,1)}
rows = []
for line in strings:
    assert line[0] == '+' and len(line) == 18
    rows.append(tuple(labels[c] for c in line[1:]))
def bracket(a,b):
    return sum(x*v + z*u for (x,z),(u,v) in zip(a,b)) % 2
def add(a,b):
    return tuple((x^u,z^v) for (x,z),(u,v) in zip(a,b))
zero = ((0,0),)*17
assert all(bracket(a,b)==0 for a in rows for b in rows)
span = {zero}
for row in rows:
    span |= {add(a,row) for a in span}
assert len(span)==4096
A = Counter(sum(p!=(0,0) for p in a) for a in span)
assert dict(A)=={0:1,8:30,10:720,12:1344,14:1680,16:321}
def syndrome(e):
    return tuple(bracket(e,row) for row in rows)
errors = [zero]
for w in (1,2):
    for support in combinations(range(17),w):
        for colors in product(((1,0),(0,1),(1,1)),repeat=w):
            e = list(zero)
            for site,color in zip(support,colors):
                e[site]=color
            errors.append(tuple(e))
assert len(errors)==1276
assert len({syndrome(e) for e in errors})==1276
# Any weight <=4 error splits into two errors of weight <=2.
# Distinct radius-two syndromes therefore prove normalizer minimum >=5.
witness = 269748239
logical = tuple(((witness>>i)&1,(witness>>(17+i))&1) for i in range(17))
assert sum(p!=(0,0) for p in logical)==5
assert syndrome(logical)==(0,)*12 and logical not in span
print('PASS: rank 12, isotropic; 1276 distinct radius-two syndromes; exact pure distance 5.')
print('Stabilizer enumerator:',dict(sorted(A.items())))
