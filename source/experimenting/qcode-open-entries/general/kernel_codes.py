"""Kernel subspaces of the prefix space E for a pure code with S_T = <p>, p = Z^m.

E = (p^perp in V_T)/<p>, dim 2m-2, coordinates: bit i (i<m-1) = Z_i, bit m-1+i = X_i X_{m-1}.
Enumerates every subspace K of E (each exactly once, via the basis b_1 < b_2 < ... with b_i the
minimum of K minus span(b_1..b_{i-1})) whose nonzero elements all have cost >= d.
Records the number of such K by (dim K, dim rad K).  Usage: python3 kernel_codes.py m [d] [maxdim]
"""
import sys, time, collections
import numpy as np

m = int(sys.argv[1]); d = int(sys.argv[2]) if len(sys.argv) > 2 else 5
maxdim = int(sys.argv[3]) if len(sys.argv) > 3 else 99
a = 2 * m - 2; N = 1 << a; ar = np.arange(N)
M = (1 << m) - 1
pop = np.vectorize(lambda v: bin(int(v)).count('1'))
basis = [1 << (m + i) for i in range(m - 1)] + [(1 << i) | (1 << (m - 1)) for i in range(m - 1)]
v = np.zeros(N, np.int64)
for j, b in enumerate(basis):
    v ^= np.where((ar >> j) & 1, b, 0)
phys = lambda u: pop((u & M) | (u >> m))
cost = np.minimum(phys(v), phys(v ^ (M << m)))
X = v & M; Z = v >> m
def symp(i, j):
    return bin(int((X[i] & Z[j]) ^ (Z[i] & X[j]))).count('1') & 1
good = cost >= d; good[0] = True
print(f'm={m} d={d}: |E|={N}, elements of cost>={d}: {int(good.sum()) - 1}', flush=True)

def rad_dim(B):
    G = [[symp(x, y) for y in B] for x in B]
    r = 0; G = [row[:] for row in G]; t = len(B)
    rows = [sum(G[i][j] << j for j in range(t)) for i in range(t)]
    piv = {}
    for x in rows:
        for b in sorted(piv, reverse=True):
            if x >> b & 1: x ^= piv[b]
        if x: piv[x.bit_length() - 1] = x
    return t - len(piv)

count = collections.Counter(); best = [0]; t0 = time.time(); nodes = [0]
def rec(B, span, A, last):
    nodes[0] += 1
    count[(len(B), rad_dim(B))] += 1
    best[0] = max(best[0], len(B))
    if len(B) == maxdim:
        return
    sp = np.array(span)
    for e in np.flatnonzero(A):
        e = int(e)
        if e <= last:
            continue
        if len(span) > 1 and (e ^ sp[1:]).min() < e:   # e must be the minimum of its coset e + span
            continue
        A2 = A & A[ar ^ e]
        rec(B + [e], span + [s ^ e for s in span], A2, e)

rec([], [0], good.copy(), 0)
print(f'max dim {best[0]}; nodes {nodes[0]}; {time.time()-t0:.0f}s')
for (t, rd), c in sorted(count.items()):
    print(f'  dim {t} rad {rd}: {c}')
