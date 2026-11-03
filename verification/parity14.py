#!/usr/bin/env python3
"""Exact modular proof: H3+H4 is odd for any binary [[14,3,d>=5]].

No solver or repository import. The signed-shadow variable T_j is an
INTEGER because it counts a genuine physical affine shadow coset.
"""
from math import comb
N = 14
SIZE = 2048
MOD = 4096
COEFF = [1,-2,-6,-17,-95,-20,-24,1108,246]

def kraw(j,w):
    return sum((-1)**q*3**(j-q)*comb(w,q)*comb(N-w,j-q)
               for q in range(max(0,j-(N-w)),min(j,w)+1))

def form_coeff(w):
    # E_0=sum H_w=2048, E_j=K_j H -2048 H_j=0 (j=1..4),
    # F_j=(-1)^w K_j H -2048 Shadow_j=0 (j=0..3).
    a=COEFF[0]+sum(COEFF[j]*kraw(j,w) for j in range(1,5))
    a+=sum(COEFF[5+j]*(-1)**w*kraw(j,w) for j in range(4))
    if 1<=w<=4:a-=SIZE*COEFF[w]
    return a

assert all(kraw(0,w)==1 for w in range(15))
assert all(kraw(j,0)==comb(14,j)*3**j for j in range(15))
weights=[form_coeff(w) for w in range(15)]
shadow=[-SIZE*COEFF[5+j] for j in range(4)]
assert all(a%MOD==0 for a in shadow),shadow
assert [i for i,a in enumerate(weights) if a%MOD] == [3,4]
assert weights[3]%MOD==weights[4]%MOD==SIZE
# E0=2048, E1..E4=0, F0..F3=0, thus the combined form
# is 2048 = 2048(H3+H4) mod 4096.
print('EXACT SHADOW-PARITY CERTIFICATE PASS')
print('Krawtchouk coefficient identities: 15/15 weights')
print('Only H3,H4 have residue 2048 modulo 4096; all affine-shadow integer coefficients vanish modulo 4096')
print('Identity: H3 + H4 is odd for every binary [[14,3,d>=5]] stabilizer')
print('Together with unconditional H3=0 and 0<=H4<=3: H4 in {1,3}')
print('All calculations use exact Python integers.')
