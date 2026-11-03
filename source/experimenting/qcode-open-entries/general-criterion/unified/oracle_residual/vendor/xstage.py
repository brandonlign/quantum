"""X-stage screen.  Given a Z-stage leaf xi (xi_0..xi_{m-2}, xi_{m-1} = sum), the X-type values
eta_k = phi(X_k X_{m-1}) satisfy the linear system  Gamma_C = sum_i (xi_i eta_i^T + eta_i xi_i^T)
(entries j<l; diagonals vanish automatically).  Its solution set is an affine space A.
  U  (unary):  for each k, some eta_k in proj_k(A) with d(eta_k + xi_A) >= d - cost(X_k X_{m-1} Z_A)
               for every A of {0..m-2}.
  P  (pairs):  arc consistency of the unary domains under the joint projections proj_{k,l}(A)
               and the conditions on X_k X_l Z_A.
Both are necessary.  Usage: python3 xstage.py n m [ids]
"""
import json, sys, time, itertools
import numpy as np
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/search')
from full import prefix
from sympl import setup
from zstage import leaves

D = 5


def solve_affine(eqs, nv):
    """eqs: list of (mask, rhs) over GF(2). Returns (particular, null basis) or None."""
    piv = {}
    for mask, r in eqs:
        for b in sorted(piv, reverse=True):
            if mask >> b & 1: mask ^= piv[b][0]; r ^= piv[b][1]
        if mask:
            b = mask.bit_length() - 1
            for q in list(piv):
                if piv[q][0] >> b & 1: piv[q] = (piv[q][0] ^ mask, piv[q][1] ^ r)
            piv[b] = (mask, r)
        elif r:
            return None
    p = 0
    for b, (mask, r) in piv.items():
        if r: p |= 1 << b
    null = []
    for fv in (j for j in range(nv) if j not in piv):
        v = 1 << fv
        for b, (mask, r) in piv.items():
            if mask >> fv & 1: v |= 1 << b
        null.append(v)
    return p, null


def proj_space(p, null, idx, b):
    """Affine projection onto the concatenated eta_k (k in idx): returns array of all values."""
    def pr(v):
        return sum(((v >> (k * b)) & ((1 << b) - 1)) << (t * b) for t, k in enumerate(idx))
    base = pr(p); gens = []
    piv = {}
    for v in null:
        x = pr(v)
        for bb in sorted(piv, reverse=True):
            if x >> bb & 1: x ^= piv[bb]
        if x: piv[x.bit_length() - 1] = x
    vals = np.array([base], np.int64)
    for g in piv.values():
        vals = np.concatenate([vals, vals ^ g])
    return vals


def screen_leaf(xs, dmin, Gam, b, m, costfn, Zbits, Xbits, pairs=True):
    xi = xs + [int(np.bitwise_xor.reduce(np.array(xs)))]
    k1 = m - 1; nv = k1 * b
    eqs = []
    for j in range(b):
        for l in range(j + 1, b):
            mask = 0
            for i in range(k1):
                if xi[i] >> j & 1: mask ^= 1 << (i * b + l)
                if xi[i] >> l & 1: mask ^= 1 << (i * b + j)
            eqs.append((mask, int(Gam[j, l])))
    sol = solve_affine(eqs, nv)
    if sol is None:
        return 'LIN'
    p, null = sol
    subs = []           # (Z-part bits of E, xi sum) for all A subset of {0..m-2}
    for A in range(1 << k1):
        s = 0
        for i in range(k1):
            if A >> i & 1: s ^= xi[i]
        subs.append((A, s))
    dom = []
    for k in range(k1):
        vals = proj_space(p, null, [k], b)
        ok = np.ones(len(vals), bool)
        for A, s in subs:
            e = Xbits[k] ^ sum(Zbits[i] for i in range(k1) if A >> i & 1)
            need = D - costfn(e)
            if need > 0: ok &= dmin[vals ^ s] >= need
        dom.append(set(int(v) for v in vals[ok]))
        if not dom[k]: return 'U'
    if not pairs: return 'PASS'
    mask_b = (1 << b) - 1
    changed = True
    while changed:
        changed = False
        for k, l in itertools.combinations(range(k1), 2):
            vals = proj_space(p, null, [k, l], b)
            ek, el = vals & mask_b, vals >> b
            ok = np.isin(ek, list(dom[k])) & np.isin(el, list(dom[l]))
            for A, s in subs:
                e = Xbits[k] ^ Xbits[l] ^ sum(Zbits[i] for i in range(k1) if A >> i & 1)
                need = D - costfn(e)
                if need > 0: ok &= dmin[ek ^ el ^ s] >= need
            nk, nl = set(int(v) for v in ek[ok]), set(int(v) for v in el[ok])
            if nk != dom[k] or nl != dom[l]:
                dom[k] &= nk; dom[l] &= nl; changed = True
                if not dom[k] or not dom[l]: return 'P'
    return 'PASS'


def run(gens, f, m, pairs=True):
    dmin, syn, Gam, b = setup(gens, f)
    basis, cost = prefix(m)
    Zbits, Xbits = basis[:m - 1], basis[m - 1:]
    tally = {}
    for xs in leaves(dmin, m):
        r = screen_leaf(xs, dmin, Gam, b, m, cost, Zbits, Xbits, pairs)
        tally[r] = tally.get(r, 0) + 1
    return tally


if __name__ == '__main__':
    n, m = int(sys.argv[1]), int(sys.argv[2])
    key = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}[m - (n - 16)]
    cells = json.load(open('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/cells.json'))[key]
    ids = range(len(cells)) if len(sys.argv) < 4 else [int(i) for i in sys.argv[3].split(',')]
    surv = []; t0 = time.time()
    for i in ids:
        t = time.time(); tal = run(cells[i], n - m, m)
        print(f'class {i}: {tal}  {time.time()-t:.1f}s', flush=True)
        if tal.get('PASS'): surv.append(i)
    print(f'n={n} m={m}: survivors {surv}  total {time.time()-t0:.0f}s')
