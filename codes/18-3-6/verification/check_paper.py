"""Check the [[18,3,6]] data printed in papers/seventeen/paper/main.tex.

Reads the generator table, the stabilizer weight enumerator, the syndrome
counts, the weight-six witness and the rational certificate of the upper bound
from the LaTeX source, and checks each one in exact arithmetic. Standard
library only; run without -O.
"""
import re
import sys
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path

if sys.flags.optimize:
    raise RuntimeError("Run without -O: this check uses assertions")

HERE = Path(__file__).resolve().parent
TEX = (HERE.parents[2] / "papers/seventeen/paper/main.tex").read_text()
N = 18
MASK = (1 << N) - 1


def enc(word):
    assert len(word) == N and set(word) <= set("IXYZ")
    v = 0
    for i, a in enumerate(word):
        k = "IXZY".index(a)
        v |= (k & 1) << i | (k >> 1) << (N + i)
    return v


def sp(a, b):
    return ((a & MASK & (b >> N)).bit_count() + ((a >> N) & b & MASK).bit_count()) % 2


def wt(a):
    return ((a & MASK) | (a >> N)).bit_count()


# Generator table: rows "\(g_i\)&\texttt{...}".
table = dict(re.findall(r"\\\(g_\{?(\d+)\}?\\\)&\\texttt\{([IXYZ]{18})\}", TEX))
assert sorted(map(int, table)) == list(range(1, 16)), sorted(table)
rows = [table[str(i)] for i in range(1, 16)]
assert (HERE / "frozen_A/generators.txt").read_text().split() == rows
G = [enc(r) for r in rows]
assert all(sp(a, b) == 0 for a in G for b in G)
S = [0]
for a in G:
    S += [b ^ a for b in S]
assert len(set(S)) == 2**15
print("PASS table: 15 independent commuting generators, equal to code A")

# Stabilizer weight enumerator (Lemma Qrank).
counts = [0] * (N + 1)
for s in S:
    counts[wt(s)] += 1
printed = dict((int(e), int(c or 1)) for c, e in
               re.findall(r"(\d*)t\^\{?(\d+)\}?", TEX[TEX.index(r"\sum_{q\in Q}"):TEX.index(r"In particular, every nonzero")]))
assert counts == [1] + [printed.get(w, 0) for w in range(1, N + 1)], (counts, printed)
print("PASS enumerator: matches the printed stabilizer weight enumerator")


# Syndromes of all errors of weight at most three.
def syndrome(e):
    return sum(sp(e, g) << i for i, g in enumerate(G))


def errors(w):
    for sites in combinations(range(N), w):
        for letters in product((1, 2, 3), repeat=w):
            e = 0
            for i, k in zip(sites, letters):
                e |= (k & 1) << i | (k >> 1) << (N + i)
            yield e


low = [syndrome(e) for w in range(3) for e in errors(w)]
assert len(low) == 1432 and len(set(low)) == 1432
three = [syndrome(e) for e in errors(3)]
assert len(three) == 22032 and not set(three) & set(low)
assert "=1432\\)" in TEX and "22{,}032" in TEX
print("PASS syndromes: 1432 distinct for weight <= 2, none repeated at weight 3")

ell = re.search(r"\\ell=\\texttt\{([IXYZ]{18})\}", TEX).group(1)
L = enc(ell)
assert wt(L) == 6 and syndrome(L) == 0 and L not in set(S)
print("PASS witness: weight six, commutes with Q, not in Q")

# Rational certificate for d_max(18,3) <= 6 (Proposition ub18).
y = list(map(F, ["-56/37", "19/1212416", "29/3637248", "5/1212416",
                 "25/14548992", "11/14548992", "1/4849664", "1/14548992"]))
for m in ["56}{37", "19}{1212416", "29}{3637248", "5}{1212416",
          "25}{14548992", "11}{14548992", "1}{4849664", "1}{14548992"]:
    assert "\\tfrac{" + m + "}" in TEX, m


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def K(j, i):
    return sum((-1)**a * 3**(j - a) * C(i, a) * C(N - i, j - a) for a in range(j + 1))


size = 2**15
# Coefficient of A_i in y0*A_0 + y1*sum(A) + sum_j y_{j+1} D_j.
lhs = [y[0] * (i == 0) + y[1] + sum(y[j + 1] * (K(j, i) - size * (i == j)) for j in range(1, 7))
       for i in range(N + 1)]
rhs = y[0] + size * y[1]
# Printed combination: low-weight terms plus multiples of the shadow sums.
low_terms = {1: 332, 2: 164, 3: 115, 4: 45, 5: 29, 6: 7, 7: 4}
shadow = [F(3, 303104), F(1, 454656), F(1, 454656), F(1, 3637248), F(1, 3637248)]
rebuilt = [F(low_terms.get(i, 0), 444) + sum(c * (-1)**i * K(j, i) for j, c in enumerate(shadow))
           for i in range(N + 1)]
assert lhs[1:] == rebuilt[1:] and lhs[0] - rebuilt[0] == 0 and rhs == -1
assert "332A_1+164A_2+115A_3+45A_4+29A_5+7A_6+4A_7" in TEX
assert all(c >= 0 for c in shadow)
print("PASS bound: printed multipliers give the printed identity equal to -1")
