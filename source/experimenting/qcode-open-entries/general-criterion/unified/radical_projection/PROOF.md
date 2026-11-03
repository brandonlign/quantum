# Necessary binary projection restriction for exact-N rank-one extensions

Assume R=<p>, p=Z^m, d>floor(m/2), and all nonzero shortened centralizer words have weight at least d. Work in E=R^perp/R with a compatible weight-admissible map M, and its symplectic adjoint tau:P->E. Fix exactly N=ker(tau), a subspace of rad(P), and replace P by P/N. Its induced alternating form is Gamma_Q. Let G=rad(P/N), of dimension g.

The induced tau is injective. For c in G, tau(c) is orthogonal to im(tau)=ker(M)^perp. Hence tau(c) lies in ker(M). Its physical lifts h and h+p are prefix-only centralizer words outside R, and both have weight at least d. At each site their weights sum to one when X(h)=0 and two when X(h)=1. Consequently

    wt(X(h)) >= 2d-m.

This X vector has even weight because h commutes with p. The map X restricted to tau(G) is injective: a nonzero pure-Z class has cost at most floor(m/2)<d and cannot lie in ker(M). Thus these X vectors constitute a binary linear [m,g,>=2d-m] code, all of whose words have even weight.

For the standard prefix basis Z_i, X_i X_(m-1), the ith X coordinate of tau(c) is the evaluation c dot xi_i. The last Z singleton image is xi_(m-1)=XOR_i<m-1 xi_i. With any basis c_1,...,c_g of G, each site therefore has tag (c_j dot xi_i)_j. For every nonzero u in F_2^g, the number of tags t with u dot t=1 must be even and at least 2d-m.

The native filter enumerates every count multiset on these 2^g tags satisfying those inequalities. It retains a partial assignment exactly when its counts are componentwise bounded by some such multiset. After adding the dependent final singleton, this is equivalent to the complete binary-code condition. This is a necessary restriction only: a retained Z assignment still needs the exact Y-stage. The filter does not assume a classified prefix code, a preferred radical basis, or a preferred order of sites. Its condition is invariant under GL(g,2), physical site permutations, and the suffix automorphisms already used for first-image reduction.

Caller requirements: the Gamma matrix is the induced form on P/N. Reduced syndrome labels are lifted to original P evaluations by dual coordinate units; restricting the original delta table to those labels gives the required minimum within N^perp/C, because every representative of such a label annihilates N. The first-image group must stabilize both C and N. `quotient_input` supplies these dual units, and `syndrome_orbits` constructs the colored physical stabilizer of N, checking linearity, distances and the induced Gram matrix on every generator. These address the two upstream conditions in the independent Opus audit.

For m=8,d=6, there are 3,10,29 admissible multisets for g=1,2,3 respectively. The implementation supports g<=3. It must be called with --gamma --avoid-radical so exact N is enforced. Interrupted runs remain LIMIT.
