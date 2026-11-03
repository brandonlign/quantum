"""Independent second implementation of the prefix-map test (paper, Section 6) as SAT (CaDiCaL).
Variables M[k][j] (k<10 syndrome bit, j<2m-2 basis vector of E).  For every coordinate vector c of
E, y_c = M c is defined by XOR chains (y_c = y_{c - topbit} xor column(topbit)).  The cost is
recomputed here from scratch (own basis construction), and for each c with t = 5 - cost(c) > 0 every
syndrome s with dmin[s] < t is forbidden by one clause.  UNSAT <=> no admissible phi."""
import os as _os
CELLS = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'census', 'cells.json')
import json, sys, time
from pysat.solvers import Cadical153
pop = lambda v: bin(v).count('1')

def own_prefix(m):
    # E basis built independently: all v in V_T commuting with Z^m, reduced mod Z^m, greedy basis
    Mm = (1 << m) - 1; P = Mm << m
    comm = lambda v: pop(v & Mm) % 2 == 0           # symp(v, Z^m) = parity of x-part
    phys = lambda v: pop((v & Mm) | (v >> m))
    piv = {P.bit_length() - 1: P}; basis = []
    for v in range(1, 1 << (2 * m)):
        if not comm(v): continue
        x = v
        for b in sorted(piv, reverse=True):
            if x >> b & 1: x ^= piv[b]
        if x: piv[x.bit_length() - 1] = x; basis.append(v)
        if len(basis) == 2 * m - 2: break
    cost = lambda v: min(phys(v), phys(v ^ P))
    return basis, cost

def own_dmin(gens, f):
    # coset minima via own syndrome map: parity-check rows = basis of {h : h.g = 0 (dot) for all g}
    n2 = 2 * f; M = (1 << f) - 1
    H = []
    piv = {}
    for h in range(1, 1 << n2):
        if any(pop(h & g) & 1 for g in gens): continue
        x = h
        for b in sorted(piv, reverse=True):
            if x >> b & 1: x ^= piv[b]
        if x: piv[x.bit_length() - 1] = x; H.append(h)
        if len(H) == n2 - len(gens): break
    syn = lambda w: sum((pop(w & h) & 1) << j for j, h in enumerate(H))
    d = [99] * (1 << len(H))
    for w in range(1 << n2):
        s = syn(w); pw = pop((w & M) | (w >> f))
        if pw < d[s]: d[s] = pw
    return d, len(H)

def sat(dmin, r, m):
    basis, cost = own_prefix(m); nb = len(basis)
    nv = [0]
    def new():
        nv[0] += 1; return nv[0]
    Mv = [[new() for j in range(nb)] for k in range(r)]
    cl = []
    y = {0: None}
    def val(c):
        v = 0
        for j in range(nb):
            if c >> j & 1: v ^= basis[j]
        return v
    for c in range(1, 1 << nb):
        hb = c.bit_length() - 1; rest = c ^ (1 << hb)
        if rest == 0: y[c] = [Mv[k][hb] for k in range(r)]; continue
        out = []
        for k in range(r):
            a, b = y[rest][k], Mv[k][hb]; z = new()
            cl += [[-z, a, b], [-z, -a, -b], [z, -a, b], [z, a, -b]]
            out.append(z)
        y[c] = out
    nlow = 0
    for c in range(1, 1 << nb):
        t = 5 - cost(val(c))
        if t <= 0: continue
        nlow += 1
        for s in range(len(dmin)):
            if dmin[s] < t:
                cl.append([(-y[c][k] if s >> k & 1 else y[c][k]) for k in range(r)])
    with Cadical153(bootstrap_with=cl) as S:
        return S.solve(), nlow, len(cl)

if __name__ == '__main__':
    cells = json.load(open(CELLS))
    KEY = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}
    for arg in sys.argv[1:]:
        m, ids = arg.split(':'); m = int(m)
        lst = range(len(cells[KEY[m]])) if ids == 'all' else [int(x) for x in ids.split(',')]
        surv = []; t0 = time.time()
        for i in lst:
            dmin, r = own_dmin(cells[KEY[m]][i], 16 - m)
            assert r == 10
            ok, nl, ncl = sat(dmin, r, m)
            print(f'm={m} class {i}: {"SAT" if ok else "UNSAT"} ({ncl} clauses, {time.time()-t0:.0f}s)', flush=True)
            if ok: surv.append(i)
        print(f'm={m}: survivors {surv}', flush=True)
