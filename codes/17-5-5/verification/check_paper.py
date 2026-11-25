"""Check the [[17,5,5]] data printed in papers/sixteen/main.tex.

Reads the generator table, the weight-five witness, the rational certificate
of the upper bound and the F4-nonlinearity remark from the LaTeX source, and
checks each one in exact arithmetic. Standard library only; run without -O.
"""
import re
import sys
from fractions import Fraction as F
from math import comb
from pathlib import Path

if sys.flags.optimize:
    raise RuntimeError("Run without -O: this check uses assertions")

HERE = Path(__file__).resolve().parent
TEX = (HERE.parents[2] / "papers/sixteen/main.tex").read_text()
N = 17
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
table = dict(re.findall(r"\\\(g_\{?(\d+)\}?\\\)&\\texttt\{([IXYZ]{17})\}", TEX))
assert sorted(map(int, table)) == list(range(1, 13)), sorted(table)
rows = [table[str(i)] for i in range(1, 13)]
assert (HERE / "generators.txt").read_text().split() == ["+" + r for r in rows]
G = [enc(r) for r in rows]
assert all(sp(a, b) == 0 for a in G for b in G)
S = [0]
for a in G:
    S += [b ^ a for b in S]
assert len(set(S)) == 4096
print("PASS table: 12 independent commuting generators, equal to generators.txt")

# Weight-five logical witness.
ell = re.search(r"\\ell=\\texttt\{([IXYZ]{17})\}", TEX).group(1)
L = enc(ell)
assert wt(L) == 5 and all(sp(L, a) == 0 for a in G) and L not in set(S)
assert min(wt(s) for s in S if s) == 8
print("PASS witness: weight five, commutes with Q, not in Q; min stabilizer weight 8")

# Rational certificate for d_max(17,5) <= 5.
y = list(map(F, ["-1001/201", "25/25728", "101/437376", "197/1749504",
                 "19/583168", "5/583168", "1/583168"]))
for m in ["1001}{201", "25}{25728", "101}{437376", "197}{1749504",
          "19}{583168", "5}{583168", "1}{583168"]:
    assert "\\tfrac{" + m + "}" in TEX, m


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def K(j, i):
    return sum((-1)**a * 3**(j - a) * C(i, a) * C(N - i, j - a) for a in range(j + 1))


cs = [y[0] * (i == 0) + y[1] + sum(y[j + 1] * (K(j, i) - 4096 * (i == j)) for j in range(1, 6))
      for i in range(N + 1)]
nums = [0, 8000, 5552, 3856, 2337, 1272, 616, 252, 81, 16, 0, 0, 1, 0, 0, 4, 9, 0]
assert cs == [F(a, 3417) for a in nums] and y[0] + 4096 * y[1] == -1
printed = "8000A_1+5552A_2+3856A_3+2337A_4+1272A_5"
assert printed in TEX and "+616A_6+252A_7+81A_8+16A_9+A_{12}+4A_{15}+9A_{16}" in TEX
print("PASS bound: printed multipliers give the printed identity equal to -1")

# Remark: Q meets V_T in exactly {0, w}.
T = {2, 3, 4, 6, 7, 9, 10, 12}
supp = lambda a: {i for i in range(N) if (a >> i) & 1 or (a >> (N + i)) & 1}
inside = [s for s in S if supp(s) <= T]
w = enc("IIZZYIXZIZZIXIIII")
assert sorted(inside) == sorted([0, w]) and "IIZZYIXZIZZIXIIII" in TEX
assert r"T=\{2,3,4,6,7,9,10,12\}" in TEX
print("PASS remark: Q meets V_T in two words, so Q is not F4-linear")
