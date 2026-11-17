# Pairwise local-unitary inequivalence of Codes A, B and C

Let Π be the orthogonal projector onto a quantum code of dimension K. For each weight w, define

A_w(Π) = K^(-2) ∑_{P: wt(P)=w} |Tr(PΠ)|²,

where P ranges over the Hermitian Pauli tensor products with one choice of phase.

For a stabilizer code, Π=|S|^(-1)∑_{g∈S}g. Pauli trace orthogonality implies that |Tr(PΠ)| is K when the binary Pauli label of P belongs to S, and zero otherwise. Hence A_w(Π) is exactly the number of stabilizer elements of weight w.

Now conjugate Π by U=U₁⊗⋯⊗Uₙ. At each coordinate, conjugation fixes the identity and acts by a real orthogonal transformation on the three-dimensional traceless Hermitian span of X,Y,Z. Therefore, for every fixed support, it orthogonally mixes the trace coefficients of the Pauli operators having precisely that support. Their sum of squares is unchanged. Summing over all supports of size w proves that A_w is invariant under arbitrary product single-qubit unitaries. Permuting qubits merely permutes supports of the same size and also preserves A_w.

The independently verified counts are:

- Code A: A₆=4.
- Code B: A₆=0 and A₇=16.
- Code C: A₆=0 and A₇=12.

Thus A is inequivalent to B and C, and B is inequivalent to C, under qubit permutations and arbitrary products of single-qubit unitaries. This is stronger than inequivalence under local Clifford transformations alone.

This does not assert inequivalence under arbitrary entangling Clifford operations. B and C were constructed from A by explicit global symplectic transvections. Such transvections are realized by Clifford conjugation, with any differences in generator signs removable by a Pauli operation.
