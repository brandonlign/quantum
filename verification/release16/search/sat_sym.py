"""SAT cross-check (sat_check.py encoding: own E basis, own coset minima, CaDiCaL) PLUS the
symmetry-breaking constraints proved sound in Lemma 10 of the paper:
  (a) xi_0 < xi_1 < ... < xi_{m-1}   (xi_i = phi(Z_i), integer order of syndrome labels)
  (b) phi(X_iX_{m-1}) < phi(X_iX_{m-1}) xor xi_i,  i < m-1
      (<=> phi(X_iX_{m-1}) has a 0 at the top set bit of xi_i).
The symmetry lemma is basis/label independent: any fixed labelling of W/C may be used, and the
order is any fixed total order on labels.  Here the labels are own_dmin's syndromes (independent of
full.syndromes), so this run also re-tests the lemma on a different labelling.
usage: sat_sym.py m:ids|m:all ... [-j workers] [-p permute_seed]"""
import os as _os
CELLS = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'census', 'cells.json')
import json, sys, time, random
from concurrent.futures import ProcessPoolExecutor
from pysat.solvers import Cadical153
from sat_check import own_prefix, own_dmin
pop = lambda v: bin(v).count('1')
KEY = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}


def coords(basis, P, v):
    """coordinates c with sum_{j in c} basis_j == v (mod P)."""
    rows = [(b, 1 << j) for j, b in enumerate(basis)] + [(P, 0)]
    piv = {}
    for x, c in rows:
        for b in sorted(piv, reverse=True):
            if x >> b & 1: x ^= piv[b][0]; c ^= piv[b][1]
        if x: piv[x.bit_length() - 1] = (x, c)
    c = 0
    for b in sorted(piv, reverse=True):
        if v >> b & 1: v ^= piv[b][0]; c ^= piv[b][1]
    assert v == 0
    return c


def build(dmin, r, m, d=5):
    basis, cost = own_prefix(m); nb = len(basis)
    P = ((1 << m) - 1) << m
    nv = [0]
    def new():
        nv[0] += 1; return nv[0]
    Mv = [[new() for j in range(nb)] for k in range(r)]
    cl = []; y = {}
    val = lambda c: __import__('functools').reduce(lambda a, j: a ^ basis[j], [j for j in range(nb) if c >> j & 1], 0)
    for c in range(1, 1 << nb):
        hb = c.bit_length() - 1; rest = c ^ (1 << hb)
        if rest == 0: y[c] = [Mv[k][hb] for k in range(r)]; continue
        out = []
        for k in range(r):
            a, b = y[rest][k], Mv[k][hb]; z = new()
            cl += [[-z, a, b], [-z, -a, -b], [z, -a, b], [z, a, -b]]; out.append(z)
        y[c] = out
    for c in range(1, 1 << nb):
        t = d - cost(val(c))
        if t <= 0: continue
        for s in range(len(dmin)):
            if dmin[s] < t:
                cl.append([(-y[c][k] if s >> k & 1 else y[c][k]) for k in range(r)])
    # symmetry breaking
    Zc = [coords(basis, P, 1 << (m + i)) for i in range(m)]
    Xc = [coords(basis, P, (1 << i) | (1 << (m - 1))) for i in range(m - 1)]
    xi = [y[c] for c in Zc]
    def lt(a, b):          # a < b as integers, bit r-1 most significant
        e = [new() for _ in range(r + 1)]   # e[k]: a,b equal on bits > k-1 ... e[r] = True
        cl.append([e[r]])
        for k in range(r - 1, -1, -1):
            # e[k] -> e[k+1] and a_k == b_k
            cl.extend([[-e[k], e[k + 1]], [-e[k], -a[k], b[k]], [-e[k], a[k], -b[k]]])
            # e[k+1] -> a_k <= b_k
            cl.append([-e[k + 1], -a[k], b[k]])
            # e[k+1] and a_k == b_k -> e[k]
            cl.extend([[-e[k + 1], a[k], b[k], e[k]], [-e[k + 1], -a[k], -b[k], e[k]]])
        cl.append([-e[0]])            # not all equal
    for i in range(m - 1): lt(xi[i], xi[i + 1])
    for i in range(m - 1):
        ph = y[Xc[i]]
        for k in range(r):
            cl.append([-xi[i][k]] + [xi[i][j] for j in range(k + 1, r)] + [-ph[k]])
    return cl


def job(a):
    m, i, seed = a
    cells = json.load(open(CELLS)); t = time.time()
    gens = cells[KEY[m]][i]
    if seed is not None:          # random coordinate permutation of C: label-independence test
        f = 16 - m; rnd = random.Random(seed * 100003 + i); pm = list(range(f)); rnd.shuffle(pm)
        gens = [sum(((g >> j & 1) << pm[j]) | ((g >> (f + j) & 1) << (f + pm[j])) for j in range(f)) for g in gens]
    dmin, r = own_dmin(gens, 16 - m); assert r == 10
    with Cadical153(bootstrap_with=build(dmin, r, m)) as S:
        ok = S.solve()
    return m, i, ok, time.time() - t


if __name__ == '__main__':
    args = sys.argv[1:]; W = 6; seed = None
    if '-j' in args: k = args.index('-j'); W = int(args[k + 1]); del args[k:k + 2]
    if '-p' in args: k = args.index('-p'); seed = int(args[k + 1]); del args[k:k + 2]
    cells = json.load(open(CELLS)); jobs = []
    for a in args:
        m, ids = a.split(':'); m = int(m)
        lst = range(len(cells[KEY[m]])) if ids == 'all' else [int(x) for x in ids.split(',')]
        jobs += [(m, i, seed) for i in lst]
    surv = []
    with ProcessPoolExecutor(W) as ex:
        for m, i, ok, dt in ex.map(job, jobs):
            print(f'm={m} class {i}: {"SAT" if ok else "UNSAT"} {dt:.0f}s', flush=True)
            if ok: surv.append((m, i))
    print('survivors', surv, flush=True)
