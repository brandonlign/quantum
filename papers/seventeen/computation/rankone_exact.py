"""Exact rank-one extension test for the [[17,3,6]] proof (R = <Z^8>).

Standard library plus numpy. Input: a 9-site additive code C of binary
dimension five (generators in the layout x | z << 9).

Setting (paper, Section 2 and Section 5): T = the 8 sites of p = Z^8,
R = <Z^8>, E = R-perp_T / R (dimension 14), W = Pauli space on the other
9 sites, P = C-perp (dimension 13), syndrome map sigma: W -> F_2^13 with
sigma(w)_j = <v_j, w> for a basis v_j of P, Gamma_jk = <v_j, v_k>, and
delta(u) = minimum weight of a word with syndrome u.

By the extension criterion (Theorem 2.1) a stabilizer with S cap V_T = R,
S-perp cap W = C and every word of S-perp outside R of weight >= d exists
if and only if there is a 13 x 14 matrix M with
    Gamma = M J M^T,  cost(e) + delta(M e) >= d  for every nonzero e in E.
This program decides that question exactly, by enumeration:

Basis of E.  z_i = Z_i + R (i = 0..6) and x_i = X_i X_7 + R (i = 0..6) form a
symplectic basis.  Every class is X_B Z_A + R with B an even subset of the 8
sites, and cost(X_B Z_A) = |B| + min(|A minus B|, |B^c minus A|).

Stage 1 (Z-half).  xi_i = M(Z_i), i = 0..7, with sum xi_i = M(Z^8) = 0.
For a proper nonempty subset A, Z_A has cost min(|A|, 8-|A|) <= 4 < d, so
    delta(sum_{i in A} xi_i) >= d - min(|A|, 8-|A|) > 0.
Site permutations of T fix R and all costs, so we may take
xi_0 < xi_1 < ... < xi_7 as integers.  Gamma = M J M^T restricted to
U = span(xi)^0 (the 6-dimensional annihilator inside the coefficient space of
P) gives c^T Gamma c' = 0 for c, c' in U: U must be Gamma-isotropic.  During
the search, after k independent images the annihilator U_k (dimension 13-k)
must contain the 6-dimensional isotropic U, which requires
    rank(Gamma restricted to U_k) <= 2 (7 - k).
When equality holds, the largest isotropic subspaces of U_k have dimension
6, so U contains the radical of Gamma on U_k and all later images lie in the
annihilator of that radical; the search restricts its candidates to it.
Optionally the search is restricted to syndromes annihilated by a given
subspace N (the kernel of tau, Section 5): im M = N^0.

Stage 2 (X-half).  With f_a dual to xi_a (f_a . xi_b = [a = b]) and u_s a
basis of U, the Gram equation is equivalent to
    eta_a := M(x_a) = base_a + Xi(Y_a),   Y_a in F_2^7,
    (Y_a)_b + (Y_b)_a = f_a^T Gamma f_b   (a != b),
where base_a is fixed by u_s . base_a = u_s^T Gamma f_a, f_b . base_a = 0.
The distance conditions for |B| = 2, 4 become: for every even B,
    s_B := sum_{a in B, a < 7} Y_a  lies in  Dom_B,
    Dom_B = { s : delta(base_B + Xi(v + s)) >= d - cost(X_B Z_supp(v)) for all v },
and |B| >= 6 gives cost >= 6 >= d (no condition).  Stage 2 is an exhaustive
backtrack over Y_0..Y_6, so a returned empty list means no extension exists
for that Z-half.
"""
import itertools
import time
import numpy as np

L = 9                      # suffix sites
FULL = (1 << 13) - 1


def wt(w):
    return ((w | (w >> L)) & 511).bit_count()


def symp(a, b):
    return (((a & 511) & (b >> L)).bit_count() + ((a >> L) & (b & 511)).bit_count()) & 1


def nullspace(rows, n):
    """Basis of {v in F_2^n : popcount(r & v) even for all rows r}."""
    piv = {}
    for r in rows:
        for p, b in piv.items():
            if r >> p & 1:
                r ^= b
        if r:
            p = r.bit_length() - 1
            for q in list(piv):
                if piv[q] >> p & 1:
                    piv[q] ^= r
            piv[p] = r
    out = []
    for f in range(n):
        if f in piv:
            continue
        v = 1 << f
        for p, b in piv.items():
            if b >> f & 1:
                v |= 1 << p
        out.append(v)
    return out


def rank(vectors):
    basis = []
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def span_words(basis):
    out = [0]
    for b in basis:
        out += [w ^ b for w in out]
    return out


def solve_dual(rows, rhs_list):
    """rows: list of 13 independent 13-bit vectors r_j.  For each rhs (13-bit
    int, bit j = required r_j . x), return the x in F_2^13."""
    n = len(rows)
    # Gaussian elimination on augmented system [rows | I]
    aug = [(rows[j], 1 << j) for j in range(n)]
    piv_rows = []
    for col in range(13):
        k = next((i for i in range(len(aug)) if aug[i][0] >> col & 1), None)
        assert k is not None, 'rows not independent'
        r = aug.pop(k)
        aug = [(a ^ r[0], b ^ r[1]) if a >> col & 1 else (a, b) for a, b in aug]
        piv_rows = [(a ^ r[0], b ^ r[1]) if a >> col & 1 else (a, b) for a, b in piv_rows]
        piv_rows.append(r)
    # now piv_rows[i] has row-part = e_{col_i} ... reduce fully
    # piv_rows are (e_col, combination) after full reduction
    inv = {}
    for a, b in piv_rows:
        assert a.bit_count() == 1
        inv[a.bit_length() - 1] = b      # x_col = sum_{j in b} rhs_j
    out = []
    for rhs in rhs_list:
        x = 0
        for col, comb in inv.items():
            if (comb & rhs).bit_count() & 1:
                x |= 1 << col
        out.append(x)
    return out


class Suffix:
    """Data of one shortened code C."""

    def __init__(self, gens, d=6):
        self.d = d
        C = span_words(gens)
        self.dim = rank(gens)
        assert self.dim == 5, self.dim
        self.dC = min(wt(w) for w in C if w)
        assert self.dC >= d, ('d(C) too small', self.dC)
        self.P = nullspace([(c >> L) | ((c & 511) << L) for c in gens], 2 * L)
        assert len(self.P) == 13
        P = self.P
        self.G = [sum(symp(P[j], P[k]) << k for k in range(13)) for j in range(13)]
        self.hull = 13 - rank(self.G)
        gens1 = []
        for i in range(L):
            for x, z in ((1, 0), (0, 1), (1, 1)):
                gens1.append(self.syn((x << i) | (z << (i + L))))
        delta = np.full(8192, -1, dtype=np.int16)
        delta[0] = 0
        fr = [0]
        while fr:
            nx = []
            for u in fr:
                for s in gens1:
                    v = u ^ s
                    if delta[v] < 0:
                        delta[v] = delta[u] + 1
                        nx.append(v)
            fr = nx
        assert (delta >= 0).all()
        self.delta = delta
        # radical of P in coefficient coordinates: c with Gamma c = 0
        self.rad = nullspace(self.G, 13)  # rows of G are symmetric: G c = 0 <=> c . G_j = 0
        assert len(self.rad) == self.hull

    def syn(self, w):
        return sum(symp(p, w) << j for j, p in enumerate(self.P))

    def gamma(self, c1, c2):
        """c1^T Gamma c2."""
        g = 0
        for j in range(13):
            if c1 >> j & 1:
                g ^= self.G[j]
        return (g & c2).bit_count() & 1


def thresholds_Z(d):
    # subset size s (1..7) of the 8 Z-images -> required delta
    return [0] + [d - min(s, 8 - s) for s in range(1, 8)]


def isotropy_state(S, Ubasis, k):
    """For U_k with basis Ubasis (dimension 13-k): return (ok, rad) where ok
    says rank(Gamma | U_k) <= 2(7-k), and rad is the radical of Gamma | U_k
    when that inequality is an equality (tight), else None.

    Tight case: the largest isotropic subspaces of U_k have dimension 6, so
    the final 6-dimensional isotropic U is one of them and contains the
    radical (U + rad is isotropic).  Hence every later image lies in rad^0.
    """
    m = len(Ubasis)
    rows = []
    for b in Ubasis:
        g = 0
        for j in range(13):
            if b >> j & 1:
                g ^= S.G[j]
        rows.append(sum(((g & b2).bit_count() & 1) << i for i, b2 in enumerate(Ubasis)))
    r = rank(rows)
    if r > 2 * (7 - k):
        return False, None
    if r < 2 * (7 - k):
        return True, None
    rad = []
    for comb in nullspace(rows, m):     # rows symmetric: kernel = radical
        v = 0
        for i in range(m):
            if comb >> i & 1:
                v ^= Ubasis[i]
        rad.append(v)
    return True, rad


def isotropic_ok(S, Ubasis, k):
    return isotropy_state(S, Ubasis, k)[0]


def restrict_annihilator(Ubasis, xi):
    """U basis -> basis of {c in U : c . xi = 0}."""
    odd = [b for b in Ubasis if (b & xi).bit_count() & 1]
    even = [b for b in Ubasis if not (b & xi).bit_count() & 1]
    if not odd:
        return None  # xi dependent on previous images (cannot happen with delta > 0)
    o = odd[0]
    return even + [b ^ o for b in odd[1:]]


def enumerate_Z(S, N=(), stats=None, forced=None, use_isotropy=True, deadline=None):
    """Yield all sorted 8-tuples (xi_0<...<xi_7) of syndromes satisfying the
    subset-sum conditions and the isotropy conditions, with every xi
    annihilated by the coefficient vectors in N.  `forced` optionally fixes a
    prefix of the tuple (used only for controls)."""
    d = S.d
    TH = thresholds_Z(d)
    delta = S.delta
    OK = {t: delta >= t for t in range(0, 8)}
    allowed = OK[TH[1]].copy()
    for n in N:
        par = np.array([(n & u).bit_count() & 1 for u in range(8192)], dtype=bool)
        allowed &= ~par
    root = np.nonzero(allowed)[0]
    if stats is not None:
        stats.setdefault('nodes', [0] * 9)
    U0 = [1 << j for j in range(13)]

    def rec(chosen, sums, cand, Ub):
        k = len(chosen)
        if stats is not None:
            stats['nodes'][k] += 1
        if deadline is not None and k <= 5 and time.time() > deadline:
            raise TimeoutError
        if k == 7:
            x7 = 0
            for v in chosen:
                x7 ^= v
            if x7 > chosen[-1] and allowed[x7]:
                yield tuple(chosen) + (x7,)
            return
        cands = cand
        if forced is not None and k < len(forced):
            cands = [c for c in cand if c == forced[k]]
        for v in cands:
            v = int(v)
            Ub2 = restrict_annihilator(Ub, v)
            if Ub2 is None:
                continue
            rad = None
            if use_isotropy:
                ok, rad = isotropy_state(S, Ub2, k + 1)
                if not ok:
                    continue
            if k == 6:
                rec_c = None
            else:
                c2 = cand[cand > v]
                # new subset sums containing v
                new = [(1, v)] + [(s + 1, t ^ v) for s, lst in sums.items() for t in lst]
                for s, t in new:
                    if len(c2) == 0:
                        break
                    c2 = c2[OK[TH[s + 1]][c2 ^ t]]
                if rad:
                    for r in rad:
                        c2 = c2[(np.bitwise_count(c2 & r) & 1) == 0]
                rec_c = c2
            sums2 = {s: list(l) for s, l in sums.items()}
            sums2.setdefault(1, []).append(v)
            for s, lst in sums.items():
                sums2.setdefault(s + 1, []).extend(t ^ v for t in lst)
            if k == 6:
                yield from rec(chosen + [v], sums2, None, Ub2)
            elif len(rec_c) or k + 1 == 7:
                yield from rec(chosen + [v], sums2, rec_c, Ub2)

    first = root
    yield from rec([], {}, first, U0)


def cost_XZ(B, A):
    """B, A: 8-bit masks; B even."""
    nb = B.bit_count()
    a = (A & ~B & 255).bit_count()
    return nb + min(a, 8 - nb - a)


def solve_X(S, xi, find_all=True, cap=None, only=None):
    """Stage 2.  xi: sorted 8-tuple.  Returns (solutions, info); each solution
    is a tuple of the 7 syndromes eta_a = M(X_a X_7).  The backtrack is
    exhaustive unless find_all is False (stop at the first solution) or `cap`
    solutions have been collected.  `only` (controls) restricts eta to a given
    7-tuple."""
    d = S.d
    delta = S.delta
    X7 = xi[:7]
    # U basis = annihilator of xi_0..xi_6
    Ub = [1 << j for j in range(13)]
    for v in X7:
        Ub = restrict_annihilator(Ub, v)
    assert len(Ub) == 6
    for a in Ub:
        for b in Ub:
            assert S.gamma(a, b) == 0
    # dual vectors f_a: f_a . xi_b = [a=b], f_a . u = 0 for a fixed complement choice:
    # choose f_a in a complement; solve with rows {xi_b} as constraints on f.
    # f is a 13-bit vector; constraints: f . xi_b = [a=b] (7 eqs), f . w_s = 0 for a
    # completion w_s of xi to a basis (6 eqs).
    comp = []
    basis = list(X7)
    for j in range(13):
        if rank(basis + [1 << j]) > len(basis):
            basis.append(1 << j)
            comp.append(1 << j)
    rows = list(X7) + comp
    f = solve_dual(rows, [1 << a for a in range(7)])
    # base_a: u_s . base = Gamma(u_s, f_a), f_b . base = 0
    rows2 = list(f) + Ub
    rhs = []
    for a in range(7):
        r = 0
        for s, u in enumerate(Ub):
            r |= S.gamma(u, f[a]) << (7 + s)
        rhs.append(r)
    base = solve_dual(rows2, rhs)
    Gm = [[S.gamma(f[a], f[b]) for b in range(7)] for a in range(7)]
    # images Xi(v) for v in F_2^7
    XiV = np.zeros(128, dtype=np.int64)
    for v in range(128):
        s = 0
        for i in range(7):
            if v >> i & 1:
                s ^= X7[i]
        XiV[v] = s
    # sanity: eta from the parametrisation satisfies Gram for random Y (checked in controls)
    # domains
    dom = {}
    for nb in (2, 4):
        for B in itertools.combinations(range(8), nb):
            Bm = sum(1 << i for i in B)
            Bp = tuple(i for i in B if i < 7)
            bB = 0
            for a in Bp:
                bB ^= base[a]
            th = np.array([d - cost_XZ(Bm, v) for v in range(128)])  # supp(v) in sites 0..6
            vals = delta[np.bitwise_xor(bB, XiV)]                     # delta(bB + Xi(w)) for w
            ok = np.zeros(128, dtype=bool)
            for s in range(128):
                ok[s] = bool((vals[np.arange(128) ^ s] >= th).all())
            dom[Bp] = ok
    sols = []
    Y = [0] * 7
    allowed_y = None
    if only is not None:
        # coordinates of the given etas in the parametrisation (None if outside it)
        inv = {int(XiV[v]): v for v in range(128)}
        allowed_y = []
        for a in range(7):
            y = inv.get(only[a] ^ base[a])
            if y is None:
                return [], dict(base=base, f=f, Ub=Ub, outside_parametrisation=a)
            allowed_y.append(y)

    def bt(a):
        if a == 7:
            sols.append(tuple(int(base[i] ^ XiV[Y[i]]) for i in range(7)))
            return (not find_all) or (cap is not None and len(sols) >= cap)
        fixed_mask = 0
        fixed_val = 0
        for b in range(a):
            fixed_mask |= 1 << b
            fixed_val |= ((Y[b] >> a & 1) ^ Gm[a][b]) << b
        for y in (range(128) if allowed_y is None else (allowed_y[a],)):
            if (y & fixed_mask) != fixed_val:
                continue
            if not dom[(a,)][y]:
                continue
            good = True
            for b in range(a):
                if not dom[(b, a)][y ^ Y[b]]:
                    good = False
                    break
            if good:
                for b, c in itertools.combinations(range(a), 2):
                    if not dom[(b, c, a)][y ^ Y[b] ^ Y[c]]:
                        good = False
                        break
            if good:
                for b, c, e in itertools.combinations(range(a), 3):
                    if not dom[(b, c, e, a)][y ^ Y[b] ^ Y[c] ^ Y[e]]:
                        good = False
                        break
            if good:
                Y[a] = y
                if bt(a + 1):
                    return True
        return False

    bt(0)
    return sols, dict(base=base, f=f, Ub=Ub, domsizes={k: int(v.sum()) for k, v in dom.items()})


def full_check(S, xi, eta):
    """Directly verify Gamma = M J M^T and cost(e) + delta(Me) >= d for all
    2^14 - 1 nonzero e, for M given by xi (8-tuple) and eta (7-tuple)."""
    d = S.d
    # Gram: Gamma(c,c') = sum_i (c.xi_i)(c'.eta_i) + (c.eta_i)(c'.xi_i)
    for j in range(13):
        for k in range(13):
            c, c2 = 1 << j, 1 << k
            g = 0
            for i in range(7):
                g ^= ((c & xi[i]).bit_count() & (c2 & eta[i]).bit_count() & 1) ^ \
                     ((c & eta[i]).bit_count() & (c2 & xi[i]).bit_count() & 1)
            if g != (S.G[j] >> k & 1):
                return False, 'gram'
    for B in range(256):
        if B.bit_count() & 1:
            continue
        hB = 0
        for a in range(7):
            if B >> a & 1:
                hB ^= eta[a]
        for A in range(256):
            if B == 0 and A in (0, 255):
                continue
            zA = 0
            for i in range(8):
                if A >> i & 1:
                    zA ^= xi[i]
            if cost_XZ(B, A) + int(S.delta[hB ^ zA]) < d:
                return False, ('distance', B, A)
    return True, None
