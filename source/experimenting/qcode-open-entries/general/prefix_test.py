"""Generic prefix-map test for a pure [[n,k,d]] code with minimum stabilizer weight m and
r = dim S_T (S_T = S ∩ V_T, T = support of a min-weight stabilizer p).

Setting (see shorten/SEARCH.md, generalised): C = {w : (0,w) in S^perp} is an additive
(f=n-m, 2^{n+k-2m+r}, >=d) code; restriction gives a linear map
    phi : E = (S_T^perp ∩ V_T)/S_T  ->  W/C ,   dim E = 2m - 2r,
with cost(v) + d(phi(v)) >= d for all v not in S_T, where cost(v) = min_{s in S_T} wt(v+s)
and d(xi) = coset minimum weight.  If no phi exists for any class of C, the case is impossible.

r=1 (S_T = <Z^m>) uses ./phisearch with the symmetry breaking proven in shorten/SEARCH.md
(the proof only needs d(xi_i) > 0 and xi_i != xi_j, i.e. d >= 3).  r=2 uses ./phisearch2
(no symmetry breaking) on an explicit S_T.

usage: prefix_test.py n k d m r cellsfile key [ids|all] [workers]"""
import sys, json, pickle, subprocess, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
pop = lambda v: bin(v).count('1')
SH = '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/shorten/'


def to_xz(word, f):                     # census layout (2 bits/site, s7_census) -> x | z<<f
    x = sum(((word >> (2 * i)) & 1) << i for i in range(f))
    z = sum(((word >> (2 * i + 1)) & 1) << i for i in range(f))
    return x | (z << f)


def coset_dmin(gens, f):
    """dmin over syndromes of W/C, W = F2^{2f}; syndrome via dot-product parity checks."""
    n2 = 2 * f; piv = {}
    for g in gens:
        x = g
        for b in sorted(piv, reverse=True):
            if x >> b & 1: x ^= piv[b]
        if x: piv[x.bit_length() - 1] = x
    # parity checks H = basis of {h : h.g = 0 for all g} (ordinary dot product)
    rows = list(piv.values()); red = {}
    for r in rows:
        for b in sorted(red, reverse=True):
            if r >> b & 1: r ^= red[b]
        if r:
            b = r.bit_length() - 1
            for bb in list(red):
                if red[bb] >> b & 1: red[bb] ^= r
            red[b] = r
    free = [j for j in range(n2) if j not in red]
    H = []
    for fc in free:
        v = 1 << fc
        for b, pr in red.items():
            if pr >> fc & 1: v |= 1 << b
        H.append(v)
    assert all(pop(h & g) % 2 == 0 for h in H for g in gens)
    col = [sum(((h >> b) & 1) << j for j, h in enumerate(H)) for b in range(n2)]
    syn = np.zeros(1 << n2, dtype=np.int32); wt = np.zeros(1 << n2, dtype=np.int8)
    M = (1 << f) - 1
    for b in range(n2):
        L = 1 << b
        syn[L:2 * L] = syn[:L] ^ col[b]
    idx = np.arange(1 << n2)
    phys = (idx & M) | (idx >> f)
    wt = np.zeros(1 << n2, dtype=np.int16)
    for b in range(f): wt += (phys >> b) & 1
    d = np.full(1 << len(H), 99, dtype=np.int16)
    np.minimum.at(d, syn, wt)
    return d.tolist(), len(H)


def prefix_input(dmin, m, r, d):
    """E basis and low list.  layout on T: x | z<<m."""
    P = ((1 << m) - 1) << m                    # Z^m
    Mm = (1 << m) - 1
    phys = lambda v: pop((v & Mm) | (v >> m))
    if r == 1:
        ST = [0, P]
        basis = [1 << (m + i) for i in range(m - 1)] + [(1 << i) | (1 << (m - 1)) for i in range(m - 1)]
        mm = m
    else:
        # S_T = <Z^m, X^m> (r=2, m even): WLOG by local Clifford (both full-support, distinct letters
        # at every site => after local Clifford p=Z^m, q=X^m).
        assert r == 2 and m % 2 == 0
        Q = Mm
        ST = [0, P, Q, P ^ Q]
        sp = lambda u, v: (pop((u & Mm) & (v >> m)) + pop((u >> m) & (v & Mm))) & 1
        # basis of {v : sp(v,P)=sp(v,Q)=0} modulo S_T
        piv = {}
        for s in (P, Q):
            x = s
            for b in sorted(piv, reverse=True):
                if x >> b & 1: x ^= piv[b]
            piv[x.bit_length() - 1] = x
        basis = []
        for v in range(1, 1 << (2 * m)):
            if sp(v, P) or sp(v, Q): continue
            x = v
            for b in sorted(piv, reverse=True):
                if x >> b & 1: x ^= piv[b]
            if x: piv[x.bit_length() - 1] = x; basis.append(v)
            if len(basis) == 2 * m - 4: break
        mm = m - 1                              # phisearch2 uses nb = 2*mm-2
    cost = lambda v: min(phys(v ^ s) for s in ST)
    nb = len(basis); low = []
    for c in range(1, 1 << nb):
        v = 0
        for j in range(nb):
            if c >> j & 1: v ^= basis[j]
        t = d - cost(v)
        if t > 0: low.append((c, t))
    return f'{mm} {len(dmin)}\n' + ' '.join(map(str, dmin)) + f'\n{len(low)}\n' + \
        '\n'.join(f'{c} {t}' for c, t in low) + '\n'


ARGS = None
def job(i):
    n, k, d, m, r, gens = ARGS[0], ARGS[1], ARGS[2], ARGS[3], ARGS[4], ARGS[5][i]
    f = n - m
    dmin, codim = coset_dmin(gens, f)
    assert codim == 2 * f - (n + k - 2 * m + r), (codim,)
    prog = SH + ('phisearch' if r == 1 else 'phisearch2')
    t = time.time()
    out = subprocess.run([prog], input=prefix_input(dmin, m, r, d), capture_output=True, text=True).stdout.strip()
    return i, out, time.time() - t


def init(a):
    global ARGS; ARGS = a


if __name__ == '__main__':
    n, k, d, m, r = map(int, sys.argv[1:6]); fn, key = sys.argv[6], sys.argv[7]
    ids = sys.argv[8] if len(sys.argv) > 8 else 'all'; W = int(sys.argv[9]) if len(sys.argv) > 9 else 6
    raw = pickle.load(open(fn, 'rb'))[key] if fn.endswith('.pkl') else json.load(open(fn))[key]
    f = n - m
    cells = [[to_xz(w, f) for w in g] for g in raw] if fn.endswith('.pkl') else raw
    lst = range(len(cells)) if ids == 'all' else [int(x) for x in ids.split(',')]
    surv = []; t0 = time.time()
    with ProcessPoolExecutor(W, initializer=init, initargs=((n, k, d, m, r, cells),)) as ex:
        for i, out, dt in ex.map(job, lst, chunksize=4):
            print(i, out, f'{dt:.1f}s', flush=True)
            if out.startswith('FOUND'): surv.append(i)
    print(f'[[{n},{k},{d}]] m={m} r={r} cell {key}: {len(lst)} classes tested, survivors {surv} ({time.time()-t0:.0f}s)', flush=True)
