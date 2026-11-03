"""Positive control for full2.search on the existing [[16,4,5]] code (pure, min stabilizer weight 6,
S_T = <p>): C has dim 9, W/C has dim 11.  The search MUST find an assignment.  We also check the
code's own phi satisfies the inequality."""
import sys
import os; os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.argv = ['x', 'q1645.txt']
exec(open('control_zcosets.py').read().split("ok = True")[0])
from full import syndromes
from full2 import search
rest = [i for i in range(n) if i not in T]
f = len(rest)
def tosuf(x, z): return sum(((x >> s) & 1) << j for j, s in enumerate(rest)) | (sum(((z >> s) & 1) << j for j, s in enumerate(rest)) << f)
# local rotation: map p's letter on each site of T to Z. Represent prefix in p-frame: for site i in T,
# letter L -> (a,b) with a = symp(L, p_i) [X-like], b = 1 if L in {p_i, other} ... use basis (q_i, p_i) with
# q_i any Pauli anticommuting with p_i: coordinates L = a q_i + b p_i.
def frame(i):
    pxi, pzi = px >> i & 1, pz >> i & 1
    q = (1, 0) if pzi else (0, 1)   # X anticommutes with Z or Y; Z anticommutes with X
    return (pxi, pzi), q
def prefcoords(x, z):
    xs = zs = 0
    for j, i in enumerate(T):
        (pxi, pzi), (qx, qz) = frame(i)
        L = (x >> i & 1, z >> i & 1)
        a = (L[0] * pzi + L[1] * pxi) % 2   # symp(L,p_i) -> q component
        # L = a q + b p
        rem = ((L[0] ^ (a * qx)), (L[1] ^ (a * qz)))
        b = 1 if rem != (0, 0) else 0
        assert rem in [(0, 0), (pxi, pzi)]
        xs |= a << j; zs |= b << j
    return xs | (zs << m)
Cg = [tosuf(x, z) for x, z in Cset]
dmin, syn, r = syndromes(Cg, f)
print('codim', r)
# own phi check
from full import prefix
basis, cost = prefix(m)
own = {}
for x, z in Sp:
    own.setdefault(prefcoords(x, z), set()).add(syn(tosuf(x & ~Tm, z & ~Tm)))
assert all(len(v) == 1 for v in own.values())
bad = [(e, s) for e, ss in own.items() for s in ss if e and e != ((1 << m) - 1) << m and cost(e) + dmin[s] < 5]
print('own phi: prefixes', len(own), 'violations', len(bad))
ok, st = search(dmin, m)
print('control search:', 'FOUND (correct)' if ok else 'NOT FOUND (BUG)', st)
