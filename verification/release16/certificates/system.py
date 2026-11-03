#!/usr/bin/env python3
"""Generic split-enumerator coset relaxation (Section 5.1 / Appendix A of the [[14,3,5]] paper),
generalised to arbitrary [[n,k,d]].  Builds integer matrices; searches a Farkas
certificate by floating LP, then rationalises and VERIFIES IT EXACTLY.
A returned 'CERT' is an exact proof of infeasibility of the relaxation;
anything else proves nothing."""
from math import comb
from fractions import Fraction
from itertools import product
import numpy as np
from scipy.optimize import linprog

def kraw(n, j, w):
    return sum((-1)**a * 3**(j-a) * comb(w, a) * comb(n-w, j-a)
               for a in range(max(0, j-(n-w)), min(j, w)+1))

def parse(p):  # Pauli string -> (x bits, z bits)
    x = sum(1 << i for i, c in enumerate(p) if c in 'XY')
    z = sum(1 << i for i, c in enumerate(p) if c in 'ZY')
    return x, z

def patterns(m, gens):
    """Coset weight histograms of R=<gens> inside R^perp on m qubits."""
    if not gens:
        return [tuple(int(i == a) for i in range(m+1)) for a in range(m+1)], None
    G = [parse(g) for g in gens]
    R = {(0, 0)}
    for g in G:
        R |= {(x ^ g[0], z ^ g[1]) for x, z in R}
    sp = lambda u, v: ((u[0] & v[1]).bit_count() + (u[1] & v[0]).bit_count()) & 1
    wt = lambda u: (u[0] | u[1]).bit_count()
    seen, pats = set(), {}
    for x in range(1 << m):
        for z in range(1 << m):
            v = (x, z)
            if v in seen or any(sp(v, r) for r in R): continue
            cos = [(x ^ r[0], z ^ r[1]) for r in R]
            seen.update(cos)
            h = [0]*(m+1)
            for c in cos: h[wt(c)] += 1
            pats[tuple(h)] = pats.get(tuple(h), 0) + 1
    return sorted(pats), pats

def build(n, k, d, blocks, weights=None, extra_cells=None):
    """blocks: list of (m_i, gens_i); suffix f = n - sum m_i.
    weights: optional dict w -> exact A_w.  Returns (Meq, beq, Gle, labels)."""
    M = 2 ** (n - k)
    ms = [m for m, _ in blocks]; f = n - sum(ms)
    pats = [patterns(m, g)[0] for m, g in blocks]
    cols = list(product(*[range(len(p)) for p in pats], range(f+1)))
    # per-block transforms of each pattern
    Tb = [[[sum(kraw(m, a, j)*P[j] for j in range(m+1)) for a in range(m+1)] for P in pl]
          for (m, _), pl in zip(blocks, pats)]
    Sb = [[[sum((-1)**j*kraw(m, a, j)*P[j] for j in range(m+1)) for a in range(m+1)] for P in pl]
          for (m, _), pl in zip(blocks, pats)]
    cells = list(product(*[range(m+1) for m in ms], range(f+1)))
    def H(c, y):
        v = int(c[-1] == y[-1])
        for i in range(len(ms)): v *= pats[i][y[i]][c[i]]
        return v
    def T(c, y):
        v = kraw(f, c[-1], y[-1])
        for i in range(len(ms)): v *= Tb[i][y[i]][c[i]]
        return v
    def SH(c, y):
        v = (-1)**y[-1] * kraw(f, c[-1], y[-1])
        for i in range(len(ms)): v *= Sb[i][y[i]][c[i]]
        return v
    eq, rhs, le, lab_eq, lab_le = [], [], [], [], []
    eq.append([H(tuple([0]*len(cells[0])), y) for y in cols]); rhs.append(1); lab_eq.append('unit')
    eq.append([sum(H(c, y) for c in cells) for y in cols]); rhs.append(M); lab_eq.append('total')
    for c in cells:
        w = sum(c)
        Tc = [T(c, y) for y in cols]; Hc = [H(c, y) for y in cols]
        if 1 <= w < d:
            eq.append([t - M*h for t, h in zip(Tc, Hc)]); rhs.append(0); lab_eq.append(('low', c))
        else:
            le.append([M*h - t for t, h in zip(Tc, Hc)]); lab_le.append(('incl', c))
            le.append([-t for t in Tc]); lab_le.append(('cnonneg', c))
        le.append([-SH(c, y) for y in cols]); lab_le.append(('shadow', c))
    for w, A in (weights or {}).items():
        eq.append([sum(H(c, y) for c in cells if sum(c) == w) for y in cols]); rhs.append(A); lab_eq.append(('weight', w))
    return np.array(eq, dtype=object), np.array(rhs, dtype=object), np.array(le, dtype=object), cols, lab_eq, lab_le

def farkas(Meq, b, G):
    """Find nu (free), lam>=0 with c = lam G + nu Meq >= 0, nu.b = -1 (float), then exact-verify."""
    ne, nl, nc = Meq.shape[0], G.shape[0], Meq.shape[1]
    Mf = Meq.astype(float); Gf = G.astype(float); bf = b.astype(float)
    sc = np.maximum(np.abs(np.vstack([Mf, Gf])).max(axis=1), 1)
    Mf /= sc[:ne, None]; bf = bf / sc[:ne]; Gf /= sc[ne:, None]
    # vars [nu (ne), lam (nl)]; constraint -(Mf^T nu + Gf^T lam) <= 0; bf.nu = -1
    A_ub = -np.hstack([Mf.T, Gf.T]); b_ub = np.zeros(nc)
    A_eq = np.hstack([bf, np.zeros(nl)])[None, :]
    bounds = [(None, None)]*ne + [(0, None)]*nl
    r = linprog(np.ones(ne+nl)*0, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[-1], bounds=bounds, method='highs')
    if r.status != 0: return None, r.message
    nu = r.x[:ne] / sc[:ne]; lam = r.x[ne:] / sc[ne:]
    for den in (1, 10**3, 10**6, 10**9, 10**12):
        nuq = [Fraction(v).limit_denominator(den) for v in nu]
        lq = [max(Fraction(0), Fraction(v).limit_denominator(den)) for v in lam]
        ok = verify(Meq, b, G, nuq, lq)
        if ok: return (nuq, lq), ok
    return None, 'rationalisation failed'

def verify(Meq, b, G, nu, lam):
    if any(l < 0 for l in lam): return False
    nc = Meq.shape[1]
    for j in range(nc):
        c = sum(nu[i]*Meq[i, j] for i in range(len(nu)) if nu[i]) + \
            sum(lam[i]*G[i, j] for i in range(len(lam)) if lam[i])
        if c < 0: return False
    val = sum(nu[i]*b[i] for i in range(len(nu)))
    return val < 0 and f'nu.rhs={val}'

def _exact_solve(A, b):
    """Solve A x = b exactly over Q (any solution, free vars = 0); None if inconsistent."""
    A = [[Fraction(v) for v in row] + [Fraction(bb)] for row, bb in zip(A, b)]
    m, n = len(A), len(A[0]) - 1
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, m) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]; A[r] = [v * inv for v in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * bb for a, bb in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == m: break
    if any(all(v == 0 for v in row[:-1]) and row[-1] != 0 for row in A): return None
    x = [Fraction(0)] * n
    for i, c in enumerate(piv): x[c] = A[i][-1]
    return x

def farkas_exact(Meq, b, G, tol=1e-9):
    """Float LP for support, then exact solve on support with tight columns forced to c=0."""
    ne, nl, nc = Meq.shape[0], G.shape[0], Meq.shape[1]
    Mf = Meq.astype(float); Gf = G.astype(float); bf = b.astype(float)
    sc = np.maximum(np.abs(np.vstack([Mf, Gf])).max(axis=1), 1)
    Mf /= sc[:ne, None]; bf = bf / sc[:ne]; Gf /= sc[ne:, None]
    A_ub = -np.hstack([Mf.T, Gf.T]); A_eq = np.hstack([bf, np.zeros(nl)])[None, :]
    bounds = [(None, None)]*ne + [(0, None)]*nl
    # minimise total multiplier mass -> sparse vertex
    cost = np.concatenate([np.zeros(ne), np.ones(nl)])
    r = linprog(cost, A_ub=A_ub, b_ub=np.zeros(nc), A_eq=A_eq, b_eq=[-1], bounds=bounds, method='highs-ds')
    if r.status != 0: return None, r.message
    x = r.x; cvec = -(A_ub @ x)
    eq_sup = [i for i in range(ne) if abs(x[i]) > tol]
    le_sup = [i for i in range(nl) if x[ne+i] > tol]
    tight = [j for j in range(nc) if abs(cvec[j]) <= 1e-9]
    # unknowns: nu_i (i in eq_sup), lam_i (i in le_sup); in original (unscaled) units
    rows, rhs = [], []
    for j in tight:
        rows.append([Meq[i, j] for i in eq_sup] + [G[i, j] for i in le_sup]); rhs.append(0)
    rows.append([b[i] for i in eq_sup] + [0]*len(le_sup)); rhs.append(-1)
    sol = _exact_solve(rows, rhs)
    if sol is None: return None, 'exact support solve inconsistent'
    nu = [Fraction(0)]*ne; lam = [Fraction(0)]*nl
    for t, i in enumerate(eq_sup): nu[i] = sol[t]
    for t, i in enumerate(le_sup): lam[i] = sol[len(eq_sup)+t]
    ok = verify(Meq, b, G, nu, lam)
    return ((nu, lam), ok) if ok else (None, 'exact candidate failed verification')
