"""Z-subset hierarchy screen for the prefix-map problem.

Setting: S pure to distance d, p = Z^m a minimum-weight stabilizer word, S_T = <p>, C the shortened
code, x_i = phi(Z_i) in W/C.  Then sum_i x_i = 0 and, for every nonempty proper A of [m],
    d(sum_{i in A} x_i) >= d - min(|A|, m - |A|).
Levels (each a necessary condition, cheapest first):
  L1  covering radius of C >= d-1
  L2  at least m distinct cosets of weight >= d-1        (pairs force the x_i distinct when d >= 3)
  L3  the full Z-subset system has a solution            (exhaustive DFS over x_0 < ... < x_{m-2})
Input is the coset-weight table dmin over F_2^r (r = dim W/C), as produced by full.syndromes.
"""
import numpy as np


def screen(dmin, m, d=5, node_limit=None):
    dmin = np.asarray(dmin); N = len(dmin); ar = np.arange(N)
    if dmin.max() < d - 1:
        return 'L1', 0
    if int((dmin >= d - 1).sum()) < m:
        return 'L2', 0
    good = {t: dmin >= t for t in range(-m, d + 1)}
    need = lambda a: d - min(a, m - a)          # required coset weight for a subset of size a
    nodes = [0]

    def rec(j, span, lo):
        # span: list of (syndrome, size) over subsets of {0..j-1}; choose x_j > lo
        nodes[0] += 1
        if node_limit and nodes[0] > node_limit:
            raise TimeoutError
        if j == m - 1:
            return True
        allowed = np.ones(N, bool); allowed[:lo + 1] = False
        for s, a in span:
            allowed &= good[need(a + 1)][s ^ ar]
            if not allowed.any():
                return False
        for x in np.flatnonzero(allowed):
            x = int(x)
            if rec(j + 1, span + [(s ^ x, a + 1) for s, a in span], x):
                return True
        return False

    return ('PASS' if rec(0, [(0, 0)], -1) else 'L3'), nodes[0]
