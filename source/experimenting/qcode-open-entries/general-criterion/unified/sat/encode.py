"""SAT encoding of one (frame, suffix class, profile) extension case for [[17,3,6]].

Unknown: the b x a map M (columns c_0..c_{2p-1} in the adapted frame basis E=Z+F).
Clauses (all are necessary conditions of a real code, so UNSAT is an exclusion):
  1. Gamma = M J M^T                                   (THEOREM.md compatibility)
  2. cost(e) + delta(M e) >= d for every nonzero e      (centralizer distance)
  3. every stabilizer word outside R has weight >= mu   (profile minimum)
  4. #{a != 0 : some word W(a,r) has weight <= 7} <= B6+B7, and
     #{a != 0 : some word W(a,r) has weight 6} <= B6   (stage3/PROOF.md budgets)
Any model is decoded and re-checked exactly; slack in (4) is handled by blocking.

Run from unified/:  python3 -B -m sat.encode FRAME CLASS P6 P7 out.cnf
"""
import json
import sys
from pathlib import Path

import numpy as np

from mixed_frames.frames import xor_table
from mixed_frames.stage3.common import (budget, graph_counts, physical_translate_costs,
                                        profile_mu, target)
from mixed_frames.stage2.common import graph_data
from model import pullback_gram


class CNF:
    def __init__(self):
        self.n = 0
        self.clauses = []

    def var(self):
        self.n += 1
        return self.n

    def add(self, c):
        self.clauses.append(c)

    native_xor = False

    def xor2(self, a, b):
        """Fresh v <-> a xor b (literals)."""
        v = self.var()
        if self.native_xor:
            self.add(('x', [-v, a, b]))
        else:
            self.add([-v, a, b]); self.add([-v, -a, -b])
            self.add([v, -a, b]); self.add([v, a, -b])
        return v

    def and2(self, a, b):
        v = self.var()
        self.add([-v, a]); self.add([-v, b]); self.add([v, -a, -b])
        return v

    def write(self, path, comments=()):
        with open(path, 'w') as fh:
            for c in comments:
                fh.write(f'c {c}\n')
            fh.write(f'p cnf {self.n} {len(self.clauses)}\n')
            for c in self.clauses:
                if isinstance(c, tuple):
                    fh.write('x' + ' '.join(map(str, c[1])) + ' 0\n')
                else:
                    fh.write(' '.join(map(str, c)) + ' 0\n')


def cubes(f):
    """Cubes (lists of (bit, value)) whose union is exactly {x : f[x]}; f is a bool array of length 2^k."""
    out = []

    def rec(g, k, prefix):
        if not g.any():
            return
        if g.all():
            out.append(list(prefix)); return
        h = len(g) // 2
        lo, hi = g[:h], g[h:]
        both = lo & hi
        bit = k - 1
        rec(both, k - 1, prefix)
        rec(lo & ~both, k - 1, prefix + [(bit, 0)])
        rec(hi & ~both, k - 1, prefix + [(bit, 1)])

    k = int(len(f)).bit_length() - 1
    rec(np.asarray(f, dtype=bool), k, [])
    return out


def forbid(cnf, bits, f, extra=None):
    """Clauses forcing vector `bits` (list of literals, bit i = bits[i]) outside set f;
    with extra literal L: (bits in f) -> L."""
    for cube in cubes(f):
        clause = [(-bits[i] if v else bits[i]) for i, v in cube]
        if extra is not None:
            clause.append(extra)
        cnf.add(clause)


def at_most(cnf, lits, k):
    """Sequential counter: at most k of lits are true."""
    if k >= len(lits):
        return
    if k == 0:
        for l in lits:
            cnf.add([-l])
        return
    s = [[cnf.var() for _ in range(k)] for _ in lits]
    for i, x in enumerate(lits):
        cnf.add([-x, s[i][0]])
        if i:
            for j in range(k):
                cnf.add([-s[i-1][j], s[i][j]])
            for j in range(1, k):
                cnf.add([-x, -s[i-1][j-1], s[i][j]])
            cnf.add([-x, -s[i-1][k-1]])


def build(name, index, profile, d=6, gamma=None, mu=None, counts=None, prefix_costs=None, native_xor=False):
    """Return (cnf, meta). Optional overrides exist only for planted-solution tests."""
    model, frame = target(name, index)
    b, p = model.b, frame.p
    a_dim = 2 * p
    Gamma = model.Gamma if gamma is None else gamma
    costs = frame.costs if prefix_costs is None else prefix_costs
    if mu is None:
        mu = profile_mu(frame, profile)
    if counts is None:
        counts = budget(frame, profile)
    B6, B7 = counts
    cnf = CNF()
    cnf.native_xor = native_xor
    # M[j][i]: bit j of column i.
    M = [[cnf.var() for _ in range(a_dim)] for _ in range(b)]
    col = lambda i: [M[j][i] for j in range(b)]
    # 1. compatibility: Gamma_jk = sum_i X_i[j]Y_i[k] + Y_i[j]X_i[k]  (X_i=col i, Y_i=col i+p)
    for j in range(b):
        for k in range(j + 1, b):
            terms = []
            for i in range(p):
                terms.append(cnf.and2(M[j][i], M[k][i + p]))
                terms.append(cnf.and2(M[j][i + p], M[k][i]))
            acc = terms[0]
            for t in terms[1:-1]:
                acc = cnf.xor2(acc, t)
            want = (Gamma[j] >> k) & 1
            last = terms[-1]
            # acc xor last == want
            if want:
                cnf.add([acc, last]); cnf.add([-acc, -last])
            else:
                cnf.add([-acc, last]); cnf.add([acc, -last])
    # 2. s_e = M e, for every e; distance condition by cost class.
    distances = np.array(model.distances)
    s = {}
    for e in range(1, 1 << a_dim):
        low = e & -e
        i = low.bit_length() - 1
        rest = e ^ low
        s[e] = col(i) if rest == 0 else [cnf.xor2(x, y) for x, y in zip(s[rest], col(i))]
    bad = {c: distances < d - c for c in set(costs[1:]) if c < d}
    for e in range(1, 1 << a_dim):
        c = costs[e]
        if c < d:
            forbid(cnf, s[e], bad[c])
    # 3/4. u_a = M^T a ; prefix class e(a) has bit (i+p)%2p = u_a[i].
    pre = physical_translate_costs(frame)
    premin = np.array([min(row) for row in pre])
    rot = np.array([sum(((u >> i) & 1) << ((i + p) % a_dim) for i in range(a_dim))
                    for u in range(1 << a_dim)])
    pmin_u = premin[rot]
    ws = graph_data(model)
    row = lambda j: M[j]
    u = {}
    light, heavy6 = [], []
    for aa in range(1, 1 << b):
        low = aa & -aa
        j = low.bit_length() - 1
        rest = aa ^ low
        u[aa] = row(j) if rest == 0 else [cnf.xor2(x, y) for x, y in zip(u[rest], row(j))]
        w = pmin_u + ws[aa]
        forbid(cnf, u[aa], w < mu)
        if (w <= 7).any() and B6 + B7 < 1 << b:
            lv = cnf.var(); light.append(lv)
            forbid(cnf, u[aa], w <= 7, extra=lv)
        if (w <= 6).any():
            hv = cnf.var(); heavy6.append(hv)
            forbid(cnf, u[aa], w <= 6, extra=hv)
    at_most(cnf, light, B6 + B7)
    at_most(cnf, heavy6, B6)
    meta = dict(frame=name, class_index=index, profile=list(profile), mu=mu, counts=[B6, B7],
                b=b, p=p, M=M, d=d)
    return cnf, meta


def decode(model_vals, meta):
    b, a_dim = meta['b'], 2 * meta['p']
    return [sum(int(model_vals[meta['M'][j][i]]) << j for j in range(b)) for i in range(a_dim)]


def check(name, index, profile, columns):
    """Independent exact check of a decoded map (adapted-frame columns)."""
    model, frame = target(name, index)
    ok_crit = model.criterion(frame.model_columns(model, columns))
    cnt = graph_counts(model, frame, columns)
    mu = profile_mu(frame, profile)
    ok_prof = cnt['minimum'] >= mu and (cnt['n6'], cnt['n7']) == tuple(budget(frame, profile))
    return dict(criterion=ok_crit, counts=cnt, profile_ok=ok_prof, valid=ok_crit and ok_prof)


if __name__ == '__main__':
    name, index, p6, p7, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    cnf, meta = build(name, index, (p6, p7), native_xor=out.endswith('.xcnf'))
    cnf.write(out, [json.dumps({k: v for k, v in meta.items() if k != 'M'})])
    Path(out + '.meta.json').write_text(json.dumps(meta))
    print(cnf.n, len(cnf.clauses))
