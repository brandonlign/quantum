# The full local kernel in the 888 frame is isotropic

Let R=<Z8,X8> be the EXACT stabilizer shortening on eight sites, in a
[[17,3,>=6]] code, and H its centralizer shortening. Every word outside R
in H has weight at least six. For every h outside R, each site sees all
four Pauli letters across h+R, so

    sum_{r in R} wt(h+r)=3*8=24.

All four summands are at least six, so all are exactly six. Thus EVERY word
of H is even: zero, the three weight-eight R words, and weight six elsewhere.
Weight parity's polar form is the Pauli symplectic form; hence H is isotropic.
Its quotient K=H/R embeds as an isotropic logical subspace, so q=dim K<=k=3.

Here a=b=12. The rank identity rank Gamma=a-2q+t, with t=0, gives
suffix hull h=2q and dim N=q. This leaves ONLY

| hull | q | dim H | dim N | Gram radical after quotient |
|---:|---:|---:|---:|---:|
|0|0|2|0|0|
|2|1|3|1|1|
|4|2|4|2|2|
|6|3|5|3|3|

For q=1, a complement h has two occurrences of each Pauli letter: the four
translate weights are 8 minus those four counts, and all are six. Site
permutation therefore gives one H class.

For q=2, use the complete checked eight-site dimension-four distance-six
prefix census. Exactly three isotropic classes contain three weight-eight
words forming a plane. That plane is necessarily R, since it contains ALL
weight-eight words. One representative of each class suffices.

For q=3, every H contains some four-dimensional subspace containing R, hence
one of those three classes after physical equivalence. Enumerate ALL cosets
in Bperp/B for each such B (256 per B), keep the ones of minimum weight six,
and add a representative. Independent site-star certificates yield ONE H
class (two extension cosets of source class 12). Every resulting H has
exactly three full words, forming R. The two known general isotropic H5
classes need not both contain a full-support plane; only one does here.

Each ordered pair of independent full R words is nonzero and distinct on
every site, so a unique single-site GL(2,2) map sends them to Z,X. Applying
these maps normalizes R to <Z8,X8>. The choice of ordered basis causes no
coverage loss: the full R plane, not an ordered basis, is the localized
stabilizer datum. Within a fixed pure-Z search, use only physical actions
that additionally preserve Z8; stronger actions are not silently used.

This is a structural reduction and local-prefix classification, NOT an
extension exclusion. It inherits the explicitly recorded complete H4 census.
No novelty claim is made. The existing 888 search remains untouched.
