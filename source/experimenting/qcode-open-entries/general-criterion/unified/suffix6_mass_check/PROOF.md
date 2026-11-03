# Dimension-six suffix census by incidence mass

Let G be the group of site permutations and independent permutations of the
three nonzero Pauli letters on nine sites. Its order is 9! 6^9. A binary
dimension-five additive code C of minimum distance at least six has one
dimension-six supercode for every nonzero coset of C whose minimum weight is
at least six. Distinct nonzero cosets give distinct supercodes over GF(2).

For each C we count these cosets exactly. Gaussian reduction gives a linear
quotient map from GF(2)^18 to GF(2)^13 with kernel C. Enumerating every physical
word of weight below six marks exactly the forbidden cosets. Unmarked labels
are all and only the allowed extensions. In particular the zero coset is
marked. The implementation checks ranks; its counts are compared against the
independent exhaustive Python physical-word coset-minimum implementation.

The number of physical incidences C < D is therefore the sum, over the complete
dimension-five census, of |G|/|Aut(C)| times the allowed coset count. Every
dimension-six binary code has exactly 63 dimension-five hyperplanes, each of
distance at least six if the code has distance at least six. Thus the total
number of physical dimension-six codes is the incidence sum divided by 63.

Independently validate each listed dimension-six representative, its physical
orbit order, and inequivalence using the site-star graph. If the sum of their
orbit sizes times 63 equals the incidence sum, there is no room for a missing
dimension-six orbit. The graph records every word and its actual nonzero Pauli
letter at every site, with distinct colors for words, letters, and site centers.
Its automorphisms are precisely physical site and letter actions preserving the
code: a physical action fixes the induced word permutation uniquely.

In addition to the total equality, enumerate all 63 hyperplanes of each of the
56 listed dimension-six representatives, identify their dimension-five classes
by independent site-star certificates, and weight each by its parent's orbit
size. For each of the 32,865 dimension-five classes separately this count must
equal its orbit size times its native allowed-extension count. This pointwise
equality prevents different class counts from cancelling in the total mass.

Completeness of dimension five is inherited from the independently replayed
suffix_mass_check/verification9.json and its explicitly recorded release-parent
dependencies. This argument strengthens the dimension-six census dependency;
it is not a quantum-code exclusion or a new assertion about the release parents.
