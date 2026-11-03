"""Full prefix-map obstruction for a pure [[16,5,5]] code whose minimum stabilizer weight is m,
with p = Z^m on T (WLOG) and S_T = <p>.
  * Restriction S^perp -> V_T has image p^perp ∩ V_T and kernel C (embedded on the other f=16-m sites);
    hence phi: (p^perp ∩ V_T) -> W/C,  v -> coset of w with (v,w) in S^perp, is well defined and
    linear, phi(p) = 0, and every (v,w) with w in phi(v) lies in S^perp.
  * Purity (every nonzero vector of S^perp has weight >= 5) gives, for v not in <p> :
        cost(v) + d(phi(v)) >= 5,   cost(v) = min(wt v, wt(v+p)),
    and d(xi) >= 5 for xi = 0 is the statement d(C) >= 5.
  * dim C = 21 - (2m-1) = 22-2m, so dim W/C = 2(16-m) - (22-2m) = 10.
Search: assign phi on the basis Z_1..Z_{m-1}, X_1X_m..X_{m-1}X_m of E=(p^perp∩V_T)/<p> by exhaustive
DFS, checking the inequality on every element of the growing span.  None found => no such code."""
import os as _os
CELLS = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'census', 'cells.json')
import itertools, json, sys, time
pop = lambda v: bin(v).count('1')

def syndromes(gens, f):
    n2 = 2 * f
    # nullspace of C under the plain dot product on 2f bits (any complement works for coset labels)
    rows = list(gens); piv = {}
    for r in rows:
        for b in sorted(piv, reverse=True):
            if r >> b & 1: r ^= piv[b]
        if r:
            b = r.bit_length() - 1
            for bb in list(piv):
                if piv[bb] >> b & 1: piv[bb] ^= r
            piv[b] = r
    free = [j for j in range(n2) if j not in piv]
    H = []
    for fcol in free:
        v = 1 << fcol
        for b, pr in piv.items():
            if pr >> fcol & 1: v |= 1 << b
        H.append(v)
    assert all(pop(h & g) % 2 == 0 for h in H for g in gens)
    M = (1 << f) - 1
    syn = lambda w: sum((pop(w & h) & 1) << j for j, h in enumerate(H))
    dmin = [99] * (1 << len(H))
    for w in range(1 << n2):
        s = syn(w); pw = pop((w & M) | (w >> f))
        if pw < dmin[s]: dmin[s] = pw
    assert dmin[0] == 0 and all(dmin[syn(g)] == 0 for g in gens)
    assert min(pw for w in range(1, 1 << n2) if syn(w) == 0 for pw in [pop((w & M) | (w >> f))]) >= 5
    return dmin, syn, len(H)

def prefix(m):
    P = ((1 << m) - 1) << m          # Z^m  (layout x | z<<m)
    Mm = (1 << m) - 1
    phys = lambda v: pop((v & Mm) | (v >> m))
    cost = lambda v: min(phys(v), phys(v ^ P))
    basis = [1 << (m + i) for i in range(m - 1)] + [(1 << i) | (1 << (m - 1)) for i in range(m - 1)]
    return basis, cost

def search(dmin, m, all_solutions=False, limit=None):
    basis, cost = prefix(m)
    nsyn = len(dmin)
    sols = []
    # span list of (e, s)
    def rec(j, span):
        if j == len(basis):
            sols.append([s for e, s in span[1:2**0]] or True); return not all_solutions
        g = basis[j]
        need_g = 5 - cost(g)
        for y in range(nsyn):
            if dmin[y] < need_g: continue
            ok = True
            for e, s in span:
                if dmin[s ^ y] + cost(e ^ g) < 5: ok = False; break
            if ok:
                if rec(j + 1, span + [(e ^ g, s ^ y) for e, s in span]): return True
        return False
    rec(0, [(0, 0)])
    return bool(sols)

if __name__ == '__main__':
    cells = json.load(open(CELLS))
    which = sys.argv[1:]  # e.g. 7:979 8:23,46  or 9:all
    for arg in which:
        m, ids = arg.split(':'); m = int(m)
        key = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}[m]
        lst = range(len(cells[key])) if ids == 'all' else [int(x) for x in ids.split(',')]
        surv = []
        t = time.time()
        for i in lst:
            dmin, syn, r = syndromes(cells[key][i], 16 - m)
            assert r == 10
            ok = search(dmin, m)
            print(f'm={m} class {i}: {"SURVIVES" if ok else "excluded"}  ({time.time()-t:.0f}s)', flush=True)
            if ok: surv.append(i)
        print(f'm={m}: survivors {surv}', flush=True)
