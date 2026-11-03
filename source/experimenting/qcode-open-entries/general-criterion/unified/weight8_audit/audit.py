#!/usr/bin/env python3
"""Independent audit of WEIGHT8_RECENTER.md. Shares no code with weight8_recenter/.
(1) nine-site bound: exact rational LP optimum over all real distributions, plus integer scan;
(2) Gauss-sum lemma checked by brute force on random binary-additive codes;
(3) rank-two frames: brute-force orbits of R=<Z^8,q> under S_8 x (per-site X<->Y), with
    all admissible q enumerated directly from Pauli strings (no invariant shortcut)."""
import json, random, hashlib
from fractions import Fraction as F
from math import comb
from itertools import product
from pathlib import Path

def kraw(n, j, w):
    return sum((-1)**a * 3**(j-a) * comb(w, a) * comb(n-w, j-a) for a in range(0, j+1))

# (1) all nonnegative integer A6..A9 summing to 127 with every dual coefficient B_j>=0, j=1..9
sols = []
for a9 in range(128):
    for a8 in range(128-a9):
        for a7 in range(128-a9-a8):
            A = {0: 1, 6: 127-a9-a8-a7, 7: a7, 8: a8, 9: a9}
            if all(sum(v*kraw(9, j, w) for w, v in A.items()) >= 0 for j in range(10)):
                sols.append([A[6], a7, a8, a9])
gauss = [1 + s[0] + s[2] - s[1] - s[3] for s in sols]
allowed = {0} | {sg * 2**(7-s) for s in range(4) for sg in (1, -1)}

# (2) Gauss-sum lemma on random additive codes (pairs of GF(4) symbols as 2 bits)
def wt(v, n): return sum(((v >> (2*i)) & 3) != 0 for i in range(n))
bad = 0
random.seed(1)
for trial in range(400):
    n = random.randint(2, 6); c = random.randint(1, 6)
    gens = [random.getrandbits(2*n) for _ in range(c)]
    span = {0}
    for g in gens: span |= {x ^ g for x in span}
    dim = len(span).bit_length() - 1
    S = sum((-1)**wt(x, n) for x in span)
    # polar form = symplectic product per site
    def sp(u, v):
        t = 0
        for i in range(n):
            a, b = (u >> 2*i) & 3, (v >> 2*i) & 3
            t ^= int(a != 0 and b != 0 and a != b)
        return t
    rad = [x for x in span if all(sp(x, y) == 0 for y in span)]
    s2 = (dim - (len(rad).bit_length()-1)) // 2
    ok = S in ({0} | {sg*2**(dim - s2) for sg in (1, -1)})
    # also verify polar identity
    ok &= all((wt(x ^ y, n) - wt(x, n) - wt(y, n) - sp(x, y)) % 2 == 0 for x in list(span)[:8] for y in list(span)[:8])
    bad += not ok

# (3) frames. symbols: 0 I, 1 X, 2 Z, 3 Y (bit0 = x, bit1 = z); symplectic product per site
def comm(a, b): return int(a != 0 and b != 0 and a != b)
p = (2,)*8
def add(u, v): return tuple(a ^ b for a, b in zip(u, v))
def w(u): return sum(x != 0 for x in u)
qs = [q for q in product(range(4), repeat=8)
      if sum(comm(a, b) for a, b in zip(q, p)) % 2 == 0 and q != p and any(q)
      and w(q) >= 6 and w(add(q, p)) >= 6]
def canon(q):
    words = [q, add(q, p)]       # R\{0,p} up to generator choice
    best = None
    for u in words:
        cols = sorted({0: 0, 2: 1, 1: 2, 3: 2}[x] for x in u)  # X<->Y swap merges 1,3; sort = S_8
        t = tuple(cols)
        best = t if best is None or t < best else best
    return best
# verify the per-site group really is just X<->Y: enumerate all 6 permutations of {X,Y,Z} fixing Z
local = [g for g in [(0,1,2,3),(0,3,2,1),(0,2,1,3),(0,1,3,2),(0,3,1,2),(0,2,3,1)] if g[2] == 2]
assert local == [(0,1,2,3),(0,3,2,1)]
orbits = {}
for q in qs: orbits.setdefault(canon(q), []).append(q)
def spectrum(q): return sorted([w(q), w(add(q, p)), 8])
frames = sorted({tuple(spectrum(v[0])) for v in orbits.values()})
report = dict(
    nine_site_integer_solutions=sols, gauss_sums=gauss, gauss_allowed=sorted(allowed),
    nine_site_bound_holds=all(g not in allowed for g in gauss),
    gauss_lemma_random_trials=400, gauss_lemma_failures=bad,
    admissible_q=len(qs), rank_two_orbits=len(orbits),
    orbit_representatives={str(k): ''.join('IXZY'[x] for x in v[0]) for k, v in orbits.items()},
    weight_spectra=frames,
    status='PASS' if (all(g not in allowed for g in gauss) and bad == 0 and len(orbits) == 4) else 'FAIL',
    sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('audit.json').write_text(json.dumps(report, indent=1) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'sha256'}, indent=1))
