# Exact full-kernel search for the 888 rank-two column

Scope: an exact eight-site stabilizer shortening R=<Z8,X8> in a hypothetical
[[17,3,6]] code. This document does not by itself exclude the target. It uses
the already proved weight-eight recentering/purity and A6+A7<=2 certificates.

## Complete local and suffix cases

PREFIX.md proves that every coset h+R outside R consists of four weight-six
words. Hence the full local centralizer H is isotropic, K=H/R has dimension
q<=3, and suffix hull h=2q, dim N=q. The complete 56-class nine-site dim6
distance-six suffix census has hulls 0,2,6 only, so q2/hull4 has no suffix.

There is one local H for q1: the four letter counts of a complementary h
are (2,2,2,2). For q3 there is one H, obtained exhaustively from the three
isotropic H4 classes containing a full-support plane by checking all 256
normalizer cosets per H4. All H have exactly three full words, forming R.
The H4 census and suffix census are explicit inherited dependencies.

prepare_suffix.py enumerates all subspaces N of rad(P) of dimension q:
three lines for hull2, 1395 three-spaces for hull6. It records complete orbits
under physically verified Aut(C). The only discarded spaces have more than
two words of weight six or seven, contradicting N subset S and A6+A7<=2.
For q1 there are 63 remaining N-orbits across 33 classes; q3 has 164 across
eight classes after screening. Hull0/q0 is left to the existing search.

## Fixed-N and fixed-H necessary Z tree

Put E=Rperp/R with its physical symplectic frame, p=6, and Q=P/N. Its Gram
Gamma has dimension 2p-q and radical q. An exact extension has an injective
isometry tau:Q -> Kperp, whose adjoint M has kernel K. Since nonzero pure-Z
classes have physical cost at most four and H outside R has weight six,
K intersects the pure-Z Lagrangian trivially. Thus the p Z images are
independent, and the evaluation of K on the Z space has rank q.

Enumerate every isomorphism rad(Gamma) -> K: 1 or 168 ordered bases for
q1/q3. Permute the Z basis (and paired complementary basis) to make its
first q radical tags independent, recomputing the COMPLETE physical cost
table in that frame. For the chosen images k_j=tau(rho_j), enforce

    rho_j dot xi_i = sp(k_j,Z_i).

Subtract the appropriate first q anchors to obtain T_i in im(Gamma) and
zprime_i in Kperp. Since K is isotropic, any Gamma-preimages c_i of T_i obey

    c_i dot T_j = sp(zprime_i,zprime_j)=0.

Different preimages differ by a radical vector and give the same scalar.
The native search enforces every such pairing and Z independence. For
EVERY nonzero subset a of the Z basis it checks

    delta(XOR_{i:a_i=1} xi_i) + min_{r in R} wt(Z_a+r) >= 6.

The `--z-costs` backend uses the exact cost indexed by the subset mask, not
an artificial site-count formula. Its header p is a dimension; internal
m=p+1 is only a recursion parameter. There is no forced singleton at the
dependent full sum. Old R1-specific xdom/projection options are forbidden.

## Full Y and exactness of N

rank2_kernel/solve.py uses the general physical frame and full uncompressed
Gamma, imposes ALL M(H)=0 equations and all Gram and distance constraints.
The exact unary kernel parities are derived in kernel_y_parity/PROOF.md.
There is no diagonal phase gauge. COMPLETE with no solution is exhaustion.

These constraints cannot silently change N. M(H)=0 implies rank(M)<=2p-q.
The full-rank radical tags imply that im(tau) has a q-dimensional radical
image; its restricted Gram has rank 2p-2q, hence rank(tau)>=2p-q. Therefore
equality holds, im(tau)=Kperp, with radical K, and ker(M)=K exactly. For
n in N, all Z evaluations vanish, so tau(n) is pure Z. Since n is in rad(P),
tau(n) is also in K. Transversality gives tau(n)=0. Dimensions give exact N.
Any FOUND result is a genuine candidate and is independently checked in
physical Pauli coordinates; no global profile filter may discard a witness.

## Outer and inner physical symmetries

Before searching Z, full physical Aut(H) preserving the R plane can reduce
radical identifications. It may mix X/Z. This is an action on COMPLETE
extensions, transporting their H,R,C,N and distances. It is used ONLY for
outer identification coverage, not as a Z lex prune. All orbit members are
recorded. For the q3 prefix, the action on K has order21, leaving eight
orbits of the 168 identifications. R is the unique full-word plane, so every
Aut(H) preserves it. The implementation additionally checks this property.

Inside a fixed identification, prefix actions must preserve Z8 (hence the
pure-Z space), preserve R and fix each k_j mod R. Small physical groups are
closed on all physical coordinate units before filtering, so products fixing
the slice are retained. Their Z actions are induced by M'(e)=M(alpha e),
exactly the matrix orientation used by the lex comparison. They preserve the
cost table, and the first-nonidentity-row compression is valid by independence.

Suffix syndrome maps S come from verified physical Aut(C,N), preserve
distances and S Gamma S^T=Gamma. The filter S^T rho_j=rho_j preserves the
identification because tau'=tau S^T. Closing the small image group before
filtering is sound even if its action has a physical kernel. Proper subgroups
would also suffice. The first threshold is d-cost(Z0), i.e. four for target
cost-two Z generators, not the R1 threshold five.

Both inner groups preserve the SAME complete-extension slice. A globally
lexicographically least Z projection in a joint orbit passes the prefix lex
tests AND has a minimum first label in its suffix orbit. Therefore combining
both reductions covers every solution. A unique joint canonical form is not
required. Outer actions that mix Z/Y are never confused with these inner tests.

## Evidence and remaining dependencies

controls.json contains genuine pure-distance-three controls for q=0..3,
including a full million-word seventeen-site physical dual check. These are
NOT distance-six target constructions. y_verification.json independently
forms the full physical-coordinate linear system for Gram, H-kernel and N
constraints, enumerates its affine space (3/1/3/6 free bits), and matches the
Y solver's complete map sets. The two smallest also match the frozen reference.
z_verification.json matches an independent generic-cost oracle on 24 artificial
cases, complete small physical trees (720 and 24 maps), and nontrivial pinned
large positives. symmetry_verification.json transports every control through
the joint physical normalization and retains a genuine completion.

The native cost-mode diff and the local averaging proof were independently
reviewed in opus_workers/rank2_backend/AUDIT.md. This does not replace complete
case coverage or the inherited prefix/suffix classifications. A production
exclusion additionally requires every H,N and identification orbit complete,
hash-pinned inputs/implementation, terminal complete statistics, and every Y
leaf completely rejected. LIMIT/ERROR/partial trees never exclude.
