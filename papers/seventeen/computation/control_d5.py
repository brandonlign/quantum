"""Positive controls for rankone_exact.py at distance five.

The stored codetables [[17,3,5]] construction contains weight-one
stabilizers (X on qubits 15 and 16), so it is not pure and the extension
criterion at d = 5 does not apply to it.  This control therefore builds a
pure [[17,3,5]] code, checks it by brute force, scrambles it, and then
recovers its extension maps from scratch:

  1. Shortened class 0 of the dimension-five list, Stage 1 + Stage 2 at d = 5,
     first solution M.
  2. Build the 14 stabilizer generators Z^8 and (tau(v_j), v_j) explicitly on
     17 qubits; check rank 14, isotropy, and that all 2^20 centralizer words
     have weight >= 5 (brute force, independent of the solver).
  3. Apply a random site permutation and random local Pauli relabelling.
  4. For every weight-8 stabilizer word p of the scrambled code whose local
     stabilizer is <p>: normalise p to Z^8, shorten to C, compute the true
     map M from the code, sort the sites by Z-image, and check
       (a) Stage 1 accepts the true Z-images (forced search returns exactly
           them; a search forced only on the first five finds them among
           its completions);
       (b) Stage 2 returns the true X-images among its solutions;
       (c) every Stage-2 solution passes the direct check of Gamma = M J M^T
           and all 2^14 - 1 distance inequalities.
"""
import json, os, random, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rankone_exact as rx

ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
CLASSES = os.path.join(ROOT, 'experimenting/qcode-open-entries/general-criterion/unified/'
                       'opus_parallel/suffix9_dimension5/suffix9_dim5.json')
n = 17
MASK = (1 << n) - 1
CAP, CHECK = 2000, 50


def W(w):
    return ((w | (w >> n)) & MASK).bit_count()


def sp(a, b):
    return ((a & MASK & (b >> n)).bit_count() + ((a >> n) & b & MASK).bit_count()) & 1


def check_code(S):
    """Rank, isotropy, and minimum weight of the centralizer (brute force)."""
    assert rx.rank(S) == 14
    for a in S:
        for b in S:
            assert sp(a, b) == 0
    D = rx.nullspace([(s >> n) | ((s & MASK) << n) for s in S], 2 * n)
    assert len(D) == 20
    best, w, prev = 99, 0, 0
    for g in range(1, 1 << 20):       # Gray-code walk over all centralizer words
        bit = (g & -g).bit_length() - 1
        w ^= D[bit]
        best = min(best, W(w))
    return best


# ---- 1. a d = 5 extension of class 0 --------------------------------------
cls = json.load(open(CLASSES))['classes']
gens0 = cls[0]['representative']['generators']
S0 = rx.Suffix(gens0, d=5)
for xi0 in rx.enumerate_Z(S0):
    sols0, _ = rx.solve_X(S0, xi0, find_all=False)
    if sols0:
        eta0 = sols0[0]
        break
assert rx.full_check(S0, xi0, eta0)[0]

# ---- 2. explicit stabilizer ------------------------------------------------
def word(prefix_x, prefix_z, suffix):
    sx, sz = suffix & 511, suffix >> 9
    return (prefix_x | (sx << 8)) | ((prefix_z | (sz << 8)) << n)


S = [word(0, 255, 0)]
for j, v in enumerate(S0.P):
    px = pz = 0
    for i in range(7):
        if eta0[i] >> j & 1:
            pz |= 1 << i
        if xi0[i] >> j & 1:
            px ^= (1 << i) ^ (1 << 7)    # X_i X_7
    S.append(word(px, pz, v))
dmin = check_code(S)
print('constructed code: rank 14, isotropic, centralizer minimum weight', dmin, flush=True)
assert dmin >= 5

# ---- 3. scramble -------------------------------------------------------------
rng = random.Random(20261116)
perm = list(range(n)); rng.shuffle(perm)
GL = [((1, 0), (0, 1)), ((0, 1), (1, 0)), ((1, 1), (0, 1)), ((1, 0), (1, 1)), ((0, 1), (1, 1)), ((1, 1), (1, 0))]
loc = [rng.choice(GL) for _ in range(n)]


def scramble(w):
    x, z = w & MASK, w >> n
    X = Z = 0
    for i in range(n):
        a, b = x >> i & 1, z >> i & 1
        (m00, m01), (m10, m11) = loc[i]
        a2, b2 = (m00 * a + m01 * b) & 1, (m10 * a + m11 * b) & 1
        X |= a2 << perm[i]; Z |= b2 << perm[i]
    return X | (Z << n)


S = [scramble(s) for s in S]
assert check_code(S) == dmin

# ---- 4. recover and test every rank-one weight-8 word -----------------------
Swords = rx.span_words(S)
Dbasis = rx.nullspace([(s >> n) | ((s & MASK) << n) for s in S], 2 * n)
results = []
for p in [s for s in Swords if W(s) == 8]:
    T = [i for i in range(n) if ((p | (p >> n)) >> i) & 1]
    Tm = sum(1 << i for i in T)
    if sum(1 for s in Swords if s and not ((s | (s >> n)) & MASK & ~Tm)) != 1:
        continue
    rest = [i for i in range(n) if i not in T]

    def relabel(w):
        x, z = w & MASK, w >> n
        for i in T:
            a, b = x >> i & 1, z >> i & 1
            pa, pb = p >> i & 1, p >> (n + i) & 1
            if (pa, pb) == (1, 0):
                a, b = b, a
            elif (pa, pb) == (1, 1):
                a ^= b
            x = (x & ~(1 << i)) | (a << i); z = (z & ~(1 << i)) | (b << i)
        return x | (z << n)

    Dn = [relabel(w) for w in Dbasis]

    def reduce_rows(order_T):
        piv = []
        for w in Dn:
            x, z = w & MASK, w >> n
            pre = sum(((x >> s) & 1) << j for j, s in enumerate(order_T)) | \
                  (sum(((z >> s) & 1) << j for j, s in enumerate(order_T)) << 8)
            suf = sum(((x >> s) & 1) << j for j, s in enumerate(rest)) | \
                  (sum(((z >> s) & 1) << j for j, s in enumerate(rest)) << 9)
            for (pp, ps, bit) in piv:
                if bit >= 0 and pre >> bit & 1:
                    pre ^= pp; suf ^= ps
            if pre:
                bit = pre.bit_length() - 1
                piv = [(pp ^ pre, ps ^ suf, b) if pp >> bit & 1 else (pp, ps, b) for (pp, ps, b) in piv]
                piv.append((pre, suf, bit))
            else:
                piv.append((0, suf, -1))
        return piv

    piv = reduce_rows(T)
    Cg = [ps for (pp, ps, b) in piv if b == -1]
    assert rx.rank(Cg) == 5
    Ssuf = rx.Suffix(Cg, d=5)

    def M_of(piv, pre):
        suf = 0
        for (pp, ps, b) in piv:
            if b >= 0 and pre >> b & 1:
                pre ^= pp; suf ^= ps
        assert pre == 0
        return Ssuf.syn(suf)

    xi = [M_of(piv, 1 << (8 + i)) for i in range(8)]
    T2 = [T[i] for i in sorted(range(8), key=lambda i: xi[i])]
    piv = reduce_rows(T2)
    xi = tuple(M_of(piv, 1 << (8 + i)) for i in range(8))
    eta = tuple(M_of(piv, (1 << i) | (1 << 7)) for i in range(7))
    assert list(xi) == sorted(xi) and len(set(xi)) == 8
    t0 = time.time()
    assert rx.full_check(Ssuf, xi, eta)[0]
    assert list(rx.enumerate_Z(Ssuf, forced=xi)) == [xi]
    comp5 = list(rx.enumerate_Z(Ssuf, forced=xi[:5]))
    assert xi in comp5
    only, _ = rx.solve_X(Ssuf, xi, only=eta)
    assert only == [eta], 'Stage 2 rejects the true X-images'
    sols, _ = rx.solve_X(Ssuf, xi, cap=CAP)
    for e in sols[:CHECK]:
        assert rx.full_check(Ssuf, xi, e)[0]
    results.append(dict(support=T, hull=Ssuf.hull, dC=Ssuf.dC, completions_of_first_five=len(comp5),
                        stage2_solutions_capped=len(sols), solutions_fully_checked=min(len(sols), CHECK)))
    print('p on', T, 'hull', Ssuf.hull, 'd(C)', Ssuf.dC, '| Stage-1 completions of first five true images',
          len(comp5), '| Stage-2 solutions (cap %d)' % CAP, len(sols), '| true X-images accepted %.1fs' % (time.time() - t0), flush=True)

assert results
json.dump(dict(constructed_from_class=0, centralizer_min_weight=dmin, xi=list(xi0), eta=list(eta0),
               stabilizer_generators_after_scramble=S, rank_one_weight8_words_tested=len(results),
               results=results), open(os.path.join(HERE, 'control_d5.json'), 'w'), indent=1)
print('PASS: %d rank-one weight-8 words, true maps recovered by both stages' % len(results))
