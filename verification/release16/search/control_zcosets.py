"""Positive control for the Z_i-coset condition.  Take an existing stabilizer code S (n=16),
find its minimum stabilizer weight m and a weight-m word p.  Locally rotate p to Z^m on T,
compute C = {w on complement : (0,w) in S^perp}, and the cosets xi_i = xi(Z_i).  For an existing
code the condition MUST hold: xi_i distinct, nonzero, sum 0, wt(Z_A)+d(xi_A) >= d."""
import sys, itertools
from collections import Counter
rows = [l.split('|') for l in open(sys.argv[1]) if l.strip()]
n = 16
g = [(sum(int(c) << i for i, c in enumerate(a.split())), sum(int(c) << i for i, c in enumerate(b.split()))) for a, b in rows]
pop = lambda v: bin(v).count('1')
def spanp(gs):
    out = [(0, 0)]
    for x, z in gs: out += [(a ^ x, b ^ z) for a, b in out]
    return out
S = spanp(g)
symp = lambda u, v: (pop(u[0] & v[1]) + pop(u[1] & v[0])) & 1
# S^perp by brute force over 2^32 is too big; compute nullspace instead
def perp(gs):
    # vectors (x,z) with symp=0 to all gens: linear system in 32 bits
    M = [(z | (x << n)) for x, z in gs]  # symp(u,(x,z)) = x_u.z + z_u.x
    # solve M * v = 0 where v=(x|z<<n): row for gen (gx,gz): coefficient on x is gz, on z is gx
    R = [(gz) | (gx << n) for gx, gz in gs]
    piv = {}; rows_ = []
    for r in R:
        for b, pr in piv.items():
            if r >> b & 1: r ^= pr
        if r:
            b = r.bit_length() - 1
            for bb in list(piv):
                if piv[bb] >> b & 1: piv[bb] ^= r
            piv[b] = r
    free = [j for j in range(2 * n) if j not in piv]
    basis = []
    for f in free:
        v = 1 << f
        for b, pr in piv.items():
            if pr >> f & 1: v |= 1 << b
        basis.append((v & ((1 << n) - 1), v >> n))
    return basis
Pb = perp(g)
print('dim S', len(g), 'dim S^perp', len(Pb))
A = Counter(pop(x | z) for x, z in S)
Sp = spanp(Pb)
B = Counter(pop(x | z) for x, z in Sp)
d = min(w for w in B if w > 0 and B[w] > A.get(w, 0))
print('stabilizer weights', sorted(A.items())[:4], 'distance', d, 'pure' if min(w for w in B if w) >= d else 'impure')
m = min(w for w in A if w)
for (px, pz) in S:
    if pop(px | pz) == m: break
T = [i for i in range(n) if (px | pz) >> i & 1]
Tm = sum(1 << i for i in T)
# local rotation sending p_i -> Z on each site i in T: we just work with "letter equal to p_i" ->Z.
# Z_i analogue at site i = Pauli on site i that commutes with p_i and is nontrivial = p_i itself.
def site(i): return ((px >> i & 1) << i, (pz >> i & 1) << i)
Cset = [v for v in Sp if (v[0] | v[1]) & Tm == 0]
print('m', m, 'T', T, '|C| = 2^%d' % (len(Cset).bit_length() - 1), 'dim S_T', (sum(1 for v in S if (v[0] | v[1]) & ~Tm == 0)).bit_length() - 1)
Cs = set(Cset)
def coset_min(i_set):
    # minimum weight of suffix w with (sum_{i in A} p_i) + w in S^perp
    vx = sum(site(i)[0] for i in i_set); vz = sum(site(i)[1] for i in i_set)
    best = 99
    for x, z in Sp:
        if ((x ^ vx) | (z ^ vz)) & Tm == 0 and (x | z) & Tm == (vx | vz) & Tm:
            best = min(best, pop((x | z) & ~Tm))
    return best
ok = True
for r in (1, 2, 3):
    for Aset in itertools.combinations(T, r):
        dm = coset_min(Aset)
        if min(r, m - r) + dm < d: ok = False; print('violated', Aset, dm)
    print('|A| =', r, 'checked')
print('control condition holds' if ok else 'CONDITION FAILS ON A REAL CODE (BUG)')
