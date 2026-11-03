"""Z-stage screen: Z-subset hierarchy + symplectic compatibility (Z-coisotropy).

For each class C, enumerate xi_0 < ... < xi_{m-2} (xi_{m-1} = sum) satisfying
    d(sum_{i in A} xi_i) >= d - min(|A|, m-|A|)   for nonempty proper A,
and keep those for which Ann(xi) is totally isotropic for the Gram matrix of C^perp.
Usage: python3 zstage.py n m [ids]
"""
import json, sys, time
import numpy as np
from sympl import setup, z_coisotropic

D = 5


def leaves(dmin, m, d=D, limit=10**7):
    N = len(dmin); ar = np.arange(N)
    good = {t: dmin >= t for t in range(-m, d + 1)}
    need = lambda a: d - min(a, m - a)
    def rec(j, span, lo, xs):
        if j == m - 1:
            yield xs; return
        allowed = np.ones(N, bool); allowed[:lo + 1] = False
        for s, a in span:
            allowed &= good[need(a + 1)][s ^ ar]
            if not allowed.any(): return
        for x in np.flatnonzero(allowed):
            x = int(x)
            yield from rec(j + 1, span + [(s ^ x, a + 1) for s, a in span], x, xs + [x])
    yield from rec(0, [(0, 0)], -1, [])


def run(gens, f, m, first=False):
    dmin, syn, Gam, b = setup(gens, f)
    nl = nc = 0
    for xs in leaves(dmin, m):
        nl += 1
        full = xs + [int(np.bitwise_xor.reduce(np.array(xs)))]
        if z_coisotropic(full, Gam, b):
            nc += 1
            if first: break
    return nl, nc


if __name__ == '__main__':
    n, m = int(sys.argv[1]), int(sys.argv[2])
    key = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}[m - (n - 16)]
    cells = json.load(open('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/cells.json'))[key]
    ids = range(len(cells)) if len(sys.argv) < 4 else [int(i) for i in sys.argv[3].split(',')]
    t0 = time.time(); surv = []
    for i in ids:
        t = time.time(); nl, nc = run(cells[i], n - m, m, first=True)
        print(f'class {i}: hierarchy leaves {nl}, coisotropic {nc}  {time.time()-t:.1f}s', flush=True)
        if nc: surv.append(i)
    print(f'n={n} m={m}: survivors {surv}  total {time.time()-t0:.0f}s')
