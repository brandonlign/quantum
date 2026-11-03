"""Independent classification of (7,2^3,>=5)_4 by coordinate form types.

Principle (proved in REPORT.md): fix a basis of a binary-dimension-3 additive
code C of length 7. Coordinate j is a linear map F2^3 -> F2^2; codeword c has a
nonzero symbol at j iff c is not annihilated by the row space V_j <= (F2^3)^*.
Local GL(2,2) on the output changes the map but not V_j, and every map with the
same V_j is reached that way. So, up to local equivalence, coordinate j is one of
15 types: V_j=0, V_j=<a> (7), V_j=<a,b> (7). Site permutations make C a multiset of
7 types. All C(21,7)=116280 multisets are enumerated; weight>=5 on all 7 nonzero
inputs and dimension 3 are imposed; survivors are realized as physical rows and
classified by the exact nauty certificate of mass_census.py (site permutation +
GL(2,2) per site + any basis change, since the graph uses all codewords).
"""
import json
from itertools import combinations_with_replacement
from pathlib import Path

OUT = Path(__file__).resolve().parent
SRC = Path('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/mass_census.py')
ns = {'__file__': str(SRC)}
exec(SRC.read_text().split('if __name__ == "__main__":')[0], ns)
graph, certificate, span, wt = ns['graph'], ns['certificate'], ns['span'], ns['wt']

types = [(0, 0)] + [(a, 0) for a in range(1, 8)] + [(a, b) for a in range(1, 8) for b in range(a + 1, 8) if b < a ^ b]
assert len(types) == 15
par = lambda f, c: (f & c).bit_count() & 1
hit = [[1 if (par(a, c) or par(b, c)) else 0 for c in range(8)] for a, b in types]

total = feasible = 0; classes = {}
for ms in combinations_with_replacement(range(15), 7):
    total += 1
    if any(sum(hit[t][c] for t in ms) < 5 for c in range(1, 8)):
        continue
    feasible += 1
    # interleaved rows: generator i, site j symbol = (bit i of a) | (bit i of b)<<1
    gens = [sum((((types[t][0] >> i) & 1) | (((types[t][1] >> i) & 1) << 1)) << (2 * j) for j, t in enumerate(ms))
            for i in range(3)]
    sp = span(gens); assert len(set(sp)) == 8 and min(wt(v, 7) for v in sp if v) >= 5
    classes.setdefault(certificate(graph(gens, 7)), gens)
assert total == 116280
mine = set(classes)
cen = json.loads((OUT / 'cells_d5_7_3.json').read_text())['7,3']['interleaved_with_aut']
theirs = {certificate(graph(g, 7)) for g, _ in cen}
res = dict(multisets=total, feasible_multisets=feasible, classes=len(mine), census_classes=len(theirs),
           certificate_sets_equal=mine == theirs)
(OUT / 'verify73_types.json').write_text(json.dumps(res, indent=2) + '\n')
print(json.dumps(res)); assert mine == theirs
