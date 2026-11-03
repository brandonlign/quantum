"""Symplectic compatibility of the prefix map.

Setting: S pure to distance d, p = Z^m, S_T = <p>, C the shortened code, P = C^perp (symplectic, in W),
b = dim P = dim W/C.  Syndromes: s(w)_j = <h_j, w> for a basis h_j of P, so W/C = F_2^b.
E has the symplectic basis Z_i, X_i X_{m-1} (i < m-1).  Write xi_i = phi(Z_i), eta_i = phi(X_i X_{m-1}).

Identity (from S = (S^perp)^perp being isotropic):
    Gamma_C = sum_i (xi_i eta_i^T + eta_i xi_i^T),        Gamma_C[j,k] = <h_j, h_k>.
Corollary (Z-coisotropy): every c in F_2^b with c.xi_i = 0 for all i < m satisfies c^T Gamma c' = 0
for every other such c'.  Equivalently D = C + <x_0, ..., x_{m-1}> is coisotropic in W.
"""
import numpy as np

pop = lambda v: bin(v).count('1')


def symp(a, b, f):
    M = (1 << f) - 1
    return (pop((a & M) & (b >> f)) + pop((a >> f) & (b & M))) & 1


def nullspace_symp(gens, f):
    """Basis of the symplectic dual of <gens> in F_2^{2f} (layout x | z << f)."""
    M = (1 << f) - 1
    sw = [((g >> f) | ((g & M) << f)) for g in gens]       # <a,w>_symp = swap(a) . w
    piv = {}
    for r in sw:
        for bb in sorted(piv, reverse=True):
            if r >> bb & 1: r ^= piv[bb]
        if r:
            bb = r.bit_length() - 1
            for q in list(piv):
                if piv[q] >> bb & 1: piv[q] ^= r
            piv[bb] = r
    out = []
    for fc in [j for j in range(2 * f) if j not in piv]:
        v = 1 << fc
        for bb, pr in piv.items():
            if pr >> fc & 1: v |= 1 << bb
        out.append(v)
    assert all(symp(h, g, f) == 0 for h in out for g in gens)
    return out


def setup(gens, f):
    """Returns dmin (coset weights by syndrome), syn, Gamma (b x b over F_2), b."""
    H = nullspace_symp(gens, f); b = len(H)
    M = (1 << f) - 1
    syn = lambda w: sum(symp(h, w, f) << j for j, h in enumerate(H))
    Gam = np.array([[symp(H[j], H[k], f) for k in range(b)] for j in range(b)], np.uint8)
    # coset weights by BFS over single-site changes
    dmin = np.full(1 << b, 99, np.int64); dmin[0] = 0
    single = [syn(((s & 1) << i) | ((s >> 1) << (f + i))) for i in range(f) for s in (1, 2, 3)]
    frontier = [0]; w = 0
    while frontier:
        w += 1; nxt = []
        for s in frontier:
            for t in single:
                u = s ^ t
                if dmin[u] > w: dmin[u] = w; nxt.append(u)
        frontier = nxt
    return dmin, syn, Gam, b


def vec(s, b):
    return np.array([(s >> j) & 1 for j in range(b)], np.uint8)


def rank2(A):
    A = A.copy() % 2; r = 0; rows, cols = A.shape
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i, c]), None)
        if p is None: continue
        A[[r, p]] = A[[p, r]]
        for i in range(rows):
            if i != r and A[i, c]: A[i] ^= A[r]
        r += 1
        if r == rows: break
    return r


def annihilator(vs, b):
    """Basis (rows) of {c in F_2^b : c . v = 0 for v in vs}."""
    if not vs: return np.eye(b, dtype=np.uint8)
    A = np.array([vec(v, b) for v in vs], np.uint8)
    # row reduce A, read off nullspace
    A = A.copy(); piv = []; r = 0
    for c in range(b):
        p = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
        if p is None: continue
        A[[r, p]] = A[[p, r]]
        for i in range(A.shape[0]):
            if i != r and A[i, c]: A[i] ^= A[r]
        piv.append(c); r += 1
    free = [c for c in range(b) if c not in piv]
    out = []
    for fc in free:
        v = np.zeros(b, np.uint8); v[fc] = 1
        for i, pc in enumerate(piv):
            if A[i, fc]: v[pc] = 1
        out.append(v)
    return np.array(out, np.uint8).reshape(len(out), b)


def z_coisotropic(xis, Gam, b):
    B = annihilator(list(xis), b)
    if len(B) == 0: return True
    return not ((B.astype(np.int64) @ Gam @ B.T.astype(np.int64)) % 2).any()
