# Stage 4: a joint physical symmetry reduction for frame 888

The extension model, census completeness assumption, exact-profile budgets,
local-plane argument and Stage-3 necessary tests remain those in
`../stage3/PROOF.md`. This proof adds a symmetry reduction, not a code exclusion.
All actions below are physical site permutations/local Cliffords, hence preserve
purity, the exact stabilizer weight histogram and the full distance criterion.

## Coefficients and syndromes are different representations

Let the independent physical columns P identify C-perp with coefficients a.
For a physical suffix automorphism g preserving C, define T by
`g(P a)=P T a`. T is invertible even when the alternating Gram form Gamma is
singular. For a physical word w its syndrome is `s(w)_i=<P_i,w>`. Therefore

` s(g w) = T^{-T} s(w) `.

Indeed `(T^T s(gw))_i=<g(P_i),g(w)>=<P_i,w>`. Set S=T^{-T}.
The suffix action on a map M=[X Y] of prefix classes to syndromes is `M -> S M`.
It preserves `Gamma=M J M^T`, since `T^T Gamma T=Gamma` implies
`S Gamma S^T=Gamma`. In particular `S Gamma=Gamma T`. No inverse Gamma is
used here; substituting a coefficient action for a syndrome action would be
incorrect in general, and agreement on im(Gamma) alone would not suffice in
the singular case.

`certificates/aut_i.json` records physical site/letter permutations, T and S.
The colored codeword/letter graph has codeword vertices and three letter vertices
forming a triangle per site; edges encode every nonzero letter of each codeword.
Its automorphisms induce site permutations and bijections of the three nonzero
Pauli symbols. Such bijections are linear symplectic maps on F2^2. Color and
incidence preserve the entire codeword set C. Even a proper subgroup would be
safe; using pynauty's generating set does not introduce a completeness assumption
needed for an exclusion. Each generator is checked for physical C and Gram
preservation, invertibility, all 4096 coefficient weights, all 4096 syndrome
coset distances, `S Gamma=Gamma T`, `T^T S=I`, and the action on a physical
representative of **every** syndrome, including those outside im(Gamma).

For each retained profile/class case, the eligible first images have
`delta(s)+2>=6`. In hull zero they additionally have
`wt(P Gamma^{-1}s)+2>=mu_out`, the unchanged Stage-3 singleton necessary test.
Both conditions are suffix invariant. Singular cases use only the distance
condition. Each case's JSON lists every eligible label once, partitioned into
orbits with minimum representatives. For every nonroot member it records a
parent and generator mapping that parent to the member. Validation checks these
edges, reachability from the minimum, disjoint union equal to eligibility, and
closure under every generator. Thus it explicitly certifies coverage, rather
than supplying an unexplained list of roots. Inverse syndrome columns allow
constructive transport from any member back to its representative. The C++ root
mask is independently compared to the certificate and frozen by hash.

## Prefix action on the affine eight-point set

R is the uniform rank-two plane on eight sites. All permutations S8 preserve R.
The uniform Clifford `X<->Y, Z->Z` also preserves R, the pure-Z Lagrangian and
all four-translate weight multisets. The seven star edges `Z_i+Z_7`, i=0..6,
are the seven nonzero circuit vectors in the quotient, with unique relation
that their sum is zero. In the actual adapted basis they have coordinates
`1,2,4,8,16,32,63`.

For injective X form the affine eight-point set

`A(X)={0,X_0,...,X_5,X_0 xor ... xor X_5}`.

Its affine span has dimension six, and its only affine dependence is the full
sum of eight points. Site permutation changes the anchor (the old site 7) and
orders the other seven vertices. These are precisely the eight translations
`A -> A+a`, a in A, followed by seven-point labelings. Conversely each such
labeling is a physical site permutation. The uniform Clifford acts trivially
on pure Z, so C2 supplies no extra X-label reduction; its Y-related diagonal
freedom is retained.

At a fixed anchor sort the seven nonzero points. Assign the smallest six to
X and the largest to the dependent XOR; any six of the circuit are independent.
Of the 7! labelings at that anchor this is unique. Then choose the smallest
sorted list over all eight anchors. There are 8!=40,320 labeled X maps per
prefix orbit: a nontrivial translation stabilizer of A would pair its points
into four copies of the same difference, giving a proper affine dependence,
contrary to the unique full dependence. (The analogous four-site fixture can
have stabilizers; it is used only to test coverage, not this factor.) C2 fixes X.

The binary enforces increasing assigned X, dependent XOR larger than the last
assigned column, and the anchor minimum. It also has a hereditary partial anchor
check. At depth j, assigned points are `{0,X_0,...,X_(j-1)}` and assigned columns
must be the smallest j final nonzero points. If translating these known points
by a known anchor already gives a sorted nonzero list smaller than the current
j-column list, adding unknown points can only lower that alternative list's
first j entries. No completion can be anchor minimal. This prune is safe;
the final check includes the dependent point and all eight anchors.

## ONE joint lexicographic-minimum argument

Do not choose unrelated representatives for the two actions. Consider the
finite collection of lists obtained from **the same** X by a suffix group
element S, then an anchor translation, then sorted labeling of A. Let L be
its lexicographic minimum. It is anchor minimal, has the largest point dependent,
and corresponds to a physical transformation of the full map, not an arbitrary
GL(6,2) change with an invented complement.

Its first nonzero point r is a suffix-orbit minimum. If some group element sent
r to t<r, apply that element to the whole affine set, then sort at its current
anchor. The resulting first point is at most t, producing a smaller list and
contradicting joint minimality. Thus BOTH implemented conditions hold for the
one joint minimum. The implementation retains all lists satisfying those
conditions; it need not remove remaining suffix ties. It never normalizes two
incompatible pieces separately and claims their intersection is safe.

A constructive witness, used in the leaf tests, is even cheaper: among all 28
pair differences of A choose one with smallest suffix-orbit minimum r. Use
its certificate path inversely to map that difference to r, then normalize
anchors. Every transformed difference is at least its own orbit minimum,
which is at least r. At an endpoint anchor, r is the smallest nonzero point;
the anchor minimum consequently also starts with r. This produces a retained
X without enumerating the full suffix group. Every pair difference is a
physical cost-two pure-Z class, so all its eligibility requirements follow
from the full Stage-3 leaf's pure-Z tests.

For a true completion apply the suffix automorphism and prefix site permutation
to the full physical stabilizer. The transformed X passes Stage 4. Complete
**all** Y in that basis using the unchanged fixed-Z solver and retain its diagonal
and complementary freedom. The uniform C2 action preserves X; all its Y maps
are included rather than removing them without a separate proof.

## Coverage of necessary Stage-3 leaves, not only true solutions

Stage-3 leaves are necessary candidates and may not have any Y completion.
Its selected unary-domain sectors (popcount <=2 in the chosen complement) are
not a prefix-invariant collection. Stage 4 disables that necessary-only filter
when canonicalization is enabled. This is a deliberate relaxation; completed
Y is still checked in full. Consequently the coverage claim is one direction:
every predecessor necessary leaf has a retained representative; Stage 4 may
also emit extra necessary leaves. It is not a claim of equality of raw lists
or of both candidate orbit sets after this relaxation.

All the remaining full-X tests are invariant. The pure-Z distance tests exhaust
the Lagrangian. Ann(X) is unchanged by a prefix basis action and transformed
contragrediently by the suffix action. Its exact pure-Z graph sector and all
four translate weights are physical and hence invariant, including singular
suffix-only coefficients. For other coefficients, full-depth caps are maxima
over an entire pure-Z fiber; a prefix action permutes these fibers and a suffix
action permutes coefficients. Thus the full cap/profile lower counts are
invariant. In the nonsingular square case all pure-Z exact counts are invariant.

Passing the full tests implies passing every partial test in any retained
labeling: cap bounds only increase when unknown bits are forgotten; partial
exact counts are subsets of the full exact sector. The full disjoint exact+
cap count dominates any cap-only guarantee on Ann(X), since its actual minimum
is no greater than the cap. The rank test is hereditary because adjoining one
annihilator dimension increases restricted alternating rank by at most two;
full Ann(X) is isotropic. Independence is preserved. Disabling unary domains
therefore makes the constructive witness a Stage-4 leaf for every Stage-3 leaf,
with the same mu/count input. Properly partial-depth diagnostic leaves do not
have this affine-set statement; coverage tests concern full X only.

## Inherited exclusions and exact completion

Only profiles (0,1),(0,2) inherit the mu7 hull-zero exclusions at indices
7,32,41,50,52. All 41 singular classes are also excluded for these profiles.
The balanced-kernel argument gives hull=2q and suffix-only dimension q. Hull two
forces a weight-six radical word. Hull six requires an eight-word suffix-only
subcode with even nonzero weights at least seven, hence distance eight; the
length-nine quaternary Plotkin bound gives size at most six. Neither conclusion
is transferred to minimum-six profiles.

Profile (0,0) is omitted entirely: saved mu8 evidence excludes all 15 hull-zero
classes and the singular obstruction excludes the other 41. `provenance.py`
replays physical tables, input/code/binary/log hashes, COMPLETE zero-leaf
summaries and disjoint legacy root coverage, writing only here. `exclusions.json`
records each omitted case and its provenance; these are inherited conclusions,
not newly completed Stage-4 searches. The retained count is
`10*2 + 56*3 = 188`, and the omitted count is `56 + 46*2 = 148`.

Post-Y acceptance is unchanged: physical compatibility/distance and minimum
outside R, then **exact** outside `(B6,B7)` equality. COMPLETE reuse requires
matching frozen computational and certificate hashes, exact context, one terminal
summary, unique leaf/completion correspondence, and recomputed physical counts
and accepted lists. LIMIT/ERROR/PENDING or diagnostic frontiers never exclude
a class. Complete coverage of all certified first roots plus exhaustive Y is
required for a new class conclusion.
