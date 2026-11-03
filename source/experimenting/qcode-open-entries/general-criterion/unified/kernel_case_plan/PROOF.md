# Residual cases and distinguished-prefix equivalence

The exact q<=3 kernel case split leaves q1/t0/hull1, q2/t2/hull1 and
q2/t0/hull3 after separating the two q3 branches. For hull1 the only N choices
are zero (q1) and its unique one-dimensional radical (q2). For hull3/q2,
enumerate every one-dimensional subspace of its three-dimensional radical:
the seven distinct nonzero words. Physical automorphisms of C act on those
seven lines; verified generator permutations and complete union-find orbits
preserve all possible exact-N cases.

The local H candidates come from the previously checked complete
kernel_capacity/profiles.json list. Equal distinguished site-star graph
certificates mean physical equivalence of PAIRS (H,p), not merely H. Coloring
the p word separately forces the physical graph isomorphism to send p to p'.
The graph's separate site and letter colors give exactly site permutations
and invertible single-site Pauli maps, as proved in suffix_mass_check.
Normalizing p and p' to Z8 therefore leaves a physical equivalence preserving
the chosen R=<Z8>. Hence one representative per distinguished certificate
is sufficient for full-H coverage, independently of radical identifications.

There are 4 q1/t0, 2 q2/t2, and **6 q2/t0** distinguished classes. The eleven
saved q2/t0 profiles contain equivalences [6,7], [9,10,11], [13,14,15].
For completeness the same routine also recovers the already-known seven
q3/t0 distinguished classes from sixteen physical p choices; q3/t2 has five.

The plan is not an exclusion. Every representative H still requires all
radical-identification orbits and a complete exact extension search. It
inherits the explicitly documented small-prefix and suffix classifications.
