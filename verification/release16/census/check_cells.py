"""Completeness of census/cells.json and emptiness of the cells used in Lemma 7 of the paper.

Runs the orbit-stabilizer mass-checked lengthening census (mass_census.py, identical to
verification/s7_mass_check.py of the [[14,3,5]] repository) twice:
  * to (10,2^10): visits (7,2^4), (8,2^6), (9,2^8), (10,2^10) and asserts the labelled mass identity
    in every cell; then checks that the classes equal the representatives in cells.json as sets of
    canonical nauty certificates;
  * to (10,2^11): the cells (10,2^11), (9,2^9), (8,2^7) are empty.
Requires pynauty==2.8.8.1.  Runtime: a few minutes."""
import json
from pathlib import Path
here = Path(__file__).resolve().parent
src = (here / 'mass_census.py').read_text()
src = src.split('if __name__ == "__main__":')[0].replace('    return cells[(maxn, target_k)]', '    return cells')
ns = {'__file__': str(here / 'mass_census.py')}; exec(src, ns)

def inter(v, n):   # cells.json layout (X bits 0..n-1, Z bits n..2n-1) -> interleaved symbols
    out = 0
    for i in range(n):
        out |= ((v >> i) & 1) << (2 * i); out |= ((v >> (n + i)) & 1) << (2 * i + 1)
    return out

cells = json.loads((here / 'cells.json').read_text())
mc = ns['census'](10, 10)
for key in ['7,4', '8,6', '9,8', '10,10']:
    n, k = map(int, key.split(','))
    mine = [[inter(v, n) for v in g] for g in cells[key]]
    for g in mine:
        sp = ns['span'](g); assert len(sp) == 1 << k and min(ns['wt'](v, n) for v in sp if v) >= 5
    A = {ns['certificate'](ns['graph'](g, n)) for g in mine}
    B = {ns['certificate'](ns['graph'](g, n)) for g, _ in mc[(n, k)]}
    print(f'cell ({n},2^{k}): {len(cells[key])} representatives, {len(A)} distinct, mass-checked census {len(B)}, equal: {A == B}', flush=True)
    assert len(A) == len(cells[key]) and A == B
mc11 = ns['census'](10, 11)
for c in [(10, 11), (9, 9), (8, 7)]:
    assert len(mc11.get(c, [])) == 0
print('cells (10,2^11), (9,2^9), (8,2^7) are empty')
print('CENSUS OK')
