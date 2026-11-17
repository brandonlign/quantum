# Two verified distance-six derivatives of Code A

Codes B and C are both pure binary [[18,3,6]] stabilizer codes. Their stabilizers have minimum weight seven, while their normalizers have minimum weight six.

Each code was checked in two complete enumerations of its 2,097,152 normalizer vectors, including a self-contained Python verifier. All 2,353,347 physical Pauli errors of weights one through five have nonzero syndrome. Every coefficient of the MacWilliams transform agrees with the enumerated normalizer.

## Exact transformations

Integers encode x | (z << 18). Pauli strings list bit-zero coordinates first. For a mask v, use the binary symplectic transvection

T_v(r) = r + ⟨r,v⟩v,

where ⟨·,·⟩ is the binary symplectic pairing.

- B is T_v(A) for v=4328808577, the Pauli word YIIIIIIYIIIIIXYIII. Its nonidentity factors occupy zero-based sites 0,7,13,14.
- C is T_v(A) for v=8607273026, the Pauli word IYIIIIYIIIIIXIIYII. Its nonidentity factors occupy zero-based sites 1,6,12,15.

Their complete binary generators are provided in each code's ROWS.txt and PAULI_GENERATORS.txt. The original Code A is unchanged.

## Verified parameters and distinguishing counts

- B: rank 15, minimum stabilizer weight 7, A₇=16, and 564 weight-six logical operators. A weight-six logical witness has integer label 34359943556.
- C: rank 15, minimum stabilizer weight 7, A₇=12, and 564 weight-six logical operators. A weight-six logical witness has integer label 34359939617.

Neither code has any normalizer element of weights one through five. Thus both have distance exactly six and are pure. Full stabilizer, normalizer and logical weight enumerators are in each STANDALONE_RESULT.json.

A, B and C are pairwise inequivalent under arbitrary products of single-qubit unitaries and qubit permutations. Their stabilizer weight enumerators distinguish them: A₆ is 4 for A and 0 for B/C, while A₇ is 16 for B and 12 for C. LOCAL_UNITARY_INEQUIVALENCE.md gives the projector-based proof. They are not claimed inequivalent under general entangling Clifford transformations.

## Reproduction

Run each code's verify_standalone.py with Python 3, without `-O`. It requires only the standard library and no other files. For a separate direct error check, compile its check_errors.cpp with a C++17 compiler and run it. To check the transformations from the included base Code A, run verify_derivatives.py from this package.

## Generator phases

The row transformation is a binary Pauli transformation; it does not by itself specify the signs of Hermitian generators. If every displayed generator of A, B and C is assigned positive phase, an exact Clifford map is P_w exp(+iπP_v/4), with v as above and:

- For B, w=17435, the Pauli word XXIXXIIIIIXIIIXIII.
- For C, w=280636, the Pauli word ZIXXXXIIIIIXIIXIII.

The Pauli correction removes the signs introduced by the bare transvection unitary. PHASE_CORRECTIONS.json and check_clifford_phases.py verify all 30 generator conjugations, using a Pauli product convention independently checked against explicit one-qubit matrices.
