# Necessary exact-profile light-count tests

All statements are conditional on the supplied suffix census, the exact
profile certificates, and the pure distance-six extension model. The 13
frame-count certificates and the integer identity `A6+A7 <= 2` are replayed
in `profile_verification.json`. No LP floating-point conclusion is used.

## Distinct physical words and exact budgets

Let `P_0,...,P_(b-1)` be the independent **physical** basis of `C^perp`.
For a compatible map `M=[X Y]`, write

`h(a)=sum a_i P_i`, `e(a)=(Y^T a, X^T a)=J M^T a`.

Every stabilizer word is

`W(a,r)=(lift_E(e(a))+r, h(a))`, with `a in F_2^b`, `r in R`.

If `W(a,r)=W(a',r')`, equality of physical suffixes gives `a=a'`, and then
equality of prefixes gives `r=r'`. This requires independence of P, not
invertibility of its Gram matrix. The `a=0` words are exactly R. All four
words for every `a!=0` are distinct, outside R; words for different a are
also distinct, even when their prefix classes coincide or vanish.

The global profile counts **words**, including nonzero words of R. Subtract
R's known counts to obtain exact outside budgets `(B6,B7)` and `L=B6+B7`:

| Frame | Global profile | B6 | B7 = L7 | L | Effective outside minimum |
|---|---|---:|---:|---:|---:|
| 688 | (1,0) | 0 | 0 | 0 | 8 |
| 688 | (1,1) | 0 | 1 | 1 | 7 |
| 688 | (2,0) | 1 | 0 | 1 | 6 |
| 888 | (0,1) | 0 | 1 | 1 | 7 |
| 888 | (0,2) | 0 | 2 | 2 | 7 |
| 888 | (1,0) | 1 | 0 | 1 | 6 |
| 888 | (1,1) | 1 | 1 | 2 | 6 |
| 888 | (2,0) | 2 | 0 | 2 | 6 |

Here 688 has R weights `(6,8,8)` and 888 has `(8,8,8)`. Profile `(0,0)` for
888 has `L=0` and effective minimum 8; the existing m8,r2 exclusion applies.
It is checked in the overlap replay rather than added to the new batch.
An exact global profile requires these exact outside counts after completion.

## The chosen plane is the entire local stabilizer

This is also necessary for the fixed shortened-code dimensions. If an extra
local stabilizer t on the same eight sites existed, its coset t+R would give
four further distinct stabilizer words outside R, with suffix zero. For 888,
all four symbols occur at every site across the R translates, so their total
weight is exactly 24. Purity forces each to have weight at least 6, hence all
four have weight 6, exceeding every allowed outside budget (at most two).

For 688, R has local rank two on six sites and rank one on two sites. Across
four translates, a rank-two site contributes total weight 3; a rank-one site
contributes at most 4. Thus the four weights sum to at most `6*3+2*4=26`.
At most one outside word may be light; even assigning that word weight 6,
their sum would be at least `6+3*8=30`. Contradiction. Therefore in both
frames `S intersect V_T=R`, including the minimum-six profiles. This uses the
exact profile bound; the older minimum-seven averaging argument alone would
not justify the minimum-six case.

## Partial-depth upper bounds

At depth j, the first j assigned X columns fix
`low_j(a)=((X_0.a),...,(X_(j-1).a))`. Define

`cap_j(low)=max c_R(z,x)` over **all** z and all x with these low bits,

`U_j(a)=cap_j(low_j(a))+wt(h(a))`.

For every completion,

`min_r wt(W(a,r)) = c_R(e(a))+wt(h(a)) <= U_j(a)`.

Thus `U<8` guarantees at least one light word for that coefficient. The
minimizing translate can depend on the completion: its existence suffices,
since different a have disjoint sets of translates. Count one guaranteed
word per such a and reject if this lower count exceeds L.

For a valid pure completion all outside weights are at least 6. Hence
`U<=6` forces the minimizing word to have weight **exactly 6**, and its
count cannot exceed B6. An upper bound of 7 alone does not prove weight 7
when weight 6 is allowed. Conversely a **lower** bound <=6 proves neither
lightness nor weight 6. These tests always use upper bounds. The separate
minimum test rejects `U<mu_out`. Coefficients omitted using
`wt(h(a))+min_low cap_j(low) >= max(mu_out,8)` cannot trigger any of these
tests; this is a safe acceleration.

Unknown translates are never individually called light. The minimum test
counts only one per a. No early equality requirement is imposed: unknown
coordinates may supply the rest of the required light words later.

## Exact translates before full Z: nonsingular square case only

If `b=2p` and Gamma is nonsingular, compatibility
`Gamma=M J M^T` makes M invertible. Therefore the unique coefficient
belonging to prefix e is

`a=Gamma^{-1} M e`.

For pure-Z `e=(z,0)` supported on the assigned j coordinates, `M e=X z` is
already known. Every R translate has the exact weight

`wt(lift_E(z,0)+r)+wt(P Gamma^{-1} X z)`.

Independent X makes the coefficients for distinct nonzero z distinct. The
DFS accumulates the counts for **newly activated** z subsets once each,
including every r, and rejects when light, weight-6, or weight-7 counts
exceed the exact case budget. This count can prune at proper depths even
when the relaxed cap lists are empty. It is never used for singular Gamma.

This exact count and the relaxed partial cap count are tested separately;
they are not added, because their coefficient sets can overlap.

## Full-Z exact sector, including singular Gamma

For `a in Ann(X)`, compatibility gives `Gamma a=X(Y^T a)`. Independent X
uniquely determines z by expressing `Gamma a` in its column span. Failure
of span membership rules out compatibility. Thus `e(a)=(z,0)` is known.
Count **every** translate `(a,r)` exactly, including `z=0` suffix-only words.

For coefficients outside Ann(X), add the one-word cap lower count only.
The implementation explicitly excludes Ann(X) coefficients from that cap
pass. These disjoint contributions may be summed safely. Reject whenever
known exact weight-6 or weight-7 counts exceed B6/B7, or the combined light
lower count exceeds L. There is no inverse-Gamma step here.

Full X does not determine all Y. Consequently it does not determine all
stabilizer weights or justify requiring global profile equality. At this
stage exact profile **upper counts** apply to the exact sector, and
necessary lower counts apply elsewhere. Missing light words cannot yet be
an exclusion.

## Completed Y

The unchanged exhaustive fixed-Z solver preserves every free diagonal and
all complementary-sector constraints, verifies compatibility, pure distance,
and the physical reconstruction. For each completed M, stage 3 independently
loops over every `a!=0` and all four r translates. Accept exactly when its
minimum outside R is sufficient and its counts are exactly `(B6,B7)`.
Adding R's known counts then gives exactly the requested global profile.
Different Y maps with the same X are all checked; singular forms and suffix-only
words are included. A limited traversal or incomplete completion never excludes.

## Frame 888 and capacity

The physical frame is `R=<0xff00,0x00ff>`. Its pure-Z quotient is a six-dimensional
Lagrangian, with nonzero costs 2 or 4, all below d=6. If a nonzero pure-Z class
mapped to syndrome zero, it would have prefix cost <6 and suffix coset distance
zero, violating the centralizer distance criterion. X is therefore injective.

For each c, distinct nonzero pure-Z classes of cost <=c require distinct
nonzero syndromes with coset distance >=6-c. This gives the capacity inequalities
in `capacity_certificate.json`. They are necessary for every hull regime; no
singular conversion from coefficients to syndromes is made.

The min7/min8 hull-zero overlap uses the separately documented balanced-kernel
argument for the full plane. It does **not** remove hull-two or hull-six classes
from any of the minimum-six cases. All 56 classes have new physical tables and
inputs for all five new 888 profiles.
