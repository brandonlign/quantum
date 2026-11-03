#!/usr/bin/env python3
"""Exact audit of the *unconditional* 14-qubit low-weight inequality.

No code-purity assumption, no imported small-code table, no solver and no
rank-specific automorphism enumeration. Rebuild each q=4 Krawtchouk coefficient
from the binomial definition and check the 15-column formal integer identity.
"""
from math import comb

N = 14
SIZE = 1 << 11

def krawtchouk(j, w):
    return sum((-1)**a * 3**(j-a) * comb(w, a) * comb(N-w, j-a)
               for a in range(max(0, j-(N-w)), min(j, w)+1))

def rhs_coefficient(w):
    delta = lambda j: int(w == j)
    return (
        -296 * (-1)**w * krawtchouk(2, w)
        -88 * (-1)**w * krawtchouk(4, w)
        -9584640 * delta(0) + 6309
        -11685888 * delta(1) -1345536 * delta(2) -4681728 * delta(3)
        +1914 * (krawtchouk(1, w)-SIZE*delta(1))
        +1318 * (krawtchouk(2, w)-SIZE*delta(2))
        +279 * (krawtchouk(3, w)-SIZE*delta(3))
        +161 * (krawtchouk(4, w)-SIZE*delta(4))
        -1081344 * delta(5) -254976 * delta(6)
        -15360 * delta(7) -44032 * delta(14)
    )

coefficients = [rhs_coefficient(w) for w in range(N + 1)]
assert coefficients == [989184 * int(w == 4) for w in range(N + 1)]

# For any [[14,3,d>=5]] isotropic stabilizer, T_j=2048 A_j (1<=j<=4),
# U_2,U_4>=0, A_w>=0, A_0=1 and sum A_w=2048.
# No A_1=A_2=A_3=0 hypothesis has been inserted.
upper = -9584640 + 6309 * SIZE
assert upper == 3336192
assert upper // 11685888 == 0     # A1=0, independent of 13q shortening
assert upper // 4681728 == 0      # A3=0, independent of 3+11 split LP
assert upper // 1345536 == 2     # A2<=2 without 2+12 split LP
budgets = [(upper - 1345536 * a2) // 989184 for a2 in range(3)]
assert budgets == [3, 2, 0]
print("PASS: 15/15 exact formal columns; no-purity low-weight inequality")
print("A1=A3=0; A2<=2; A4<=3,2,0 for A2=0,1,2 respectively")
