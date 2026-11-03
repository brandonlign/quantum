"""Vectorized version of full.py (same mathematics; see full.py docstring), with a positive control.
phi is assigned on the basis of E = (p^perp ∩ V_T)/<p>, one basis vector at a time; for the new
basis vector g the allowed syndromes y satisfy dmin[s ^ y] >= 5 - cost(e ^ g) for all (e,s) in span."""
import os as _os
CELLS = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'census', 'cells.json')
import json, sys, time
import numpy as np
from full import syndromes, prefix

def search(dmin, m, d=5, count=False, node_limit=None):
    dmin = np.array(dmin); N = len(dmin); ar = np.arange(N)
    basis, cost = prefix(m)
    good = {t: dmin >= t for t in range(-6, 12)}
    stats = dict(nodes=0, sols=0)
    def rec(j, E, Sy):
        stats['nodes'] += 1
        if node_limit and stats['nodes'] > node_limit: raise TimeoutError
        if j == len(basis):
            stats['sols'] += 1; return not count
        g = basis[j]
        allowed = np.ones(N, bool)
        for e, s in zip(E, Sy):
            allowed &= good[d - cost(e ^ g)][s ^ ar]
            if not allowed.any(): return False
        for y in np.flatnonzero(allowed):
            if rec(j + 1, E + [e ^ g for e in E], Sy + [s ^ int(y) for s in Sy]): return True
        return False
    found = rec(0, [0], [0])
    return found, stats

if __name__ == '__main__':
    cells = json.load(open(CELLS))
    for arg in sys.argv[1:]:
        m, ids = arg.split(':'); m = int(m)
        key = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}[m]
        lst = range(len(cells[key])) if ids == 'all' else [int(x) for x in ids.split(',')]
        surv = []; t = time.time()
        for i in lst:
            dmin, syn, r = syndromes(cells[key][i], 16 - m)
            assert r == 10
            ok, st = search(dmin, m)
            print(f'm={m} class {i}: {"SURVIVES" if ok else "excluded"} nodes={st["nodes"]} ({time.time()-t:.0f}s)', flush=True)
            if ok: surv.append(i)
        print(f'm={m}: survivors {surv}', flush=True)
