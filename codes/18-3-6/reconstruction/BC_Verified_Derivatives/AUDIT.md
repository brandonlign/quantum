# Verified B/C-style derivatives: independent PASS

The sealed Code A was copied from VERIFIED_CODE.json without modification. The initially requested Z transvection on zero-based sites2,3,4,5 was fully enumerated and rejected: it gives minS7 and A7=16, but twelve weight-five logical operators, hence distance five.

A complete finite scan of 3,060 weight-four Z masks found no derivative having both minS7 and distance at least six. A subsequent bounded scan of general weight-four Pauli masks stopped after 86,085 of the 247,860 possible masks, as soon as both requested fingerprints were obtained. It found two valid derivatives in that prefix. This is not an exhaustive classification of transvections or codes.

The screening source uses the elementary bound |wt(x+v)-wt(x)|≤4: only original stabilizers of weight≤11 can affect the new minimum below seven or A7; only original normalizer words of weight≤10 can affect minimum below six or D6. Its frozen base normalizer is independently computed. More importantly, neither positive conclusion relies on screening: each selected derivative received a fresh nullspace, a full C++ normalizer enumeration, a self-contained second full Python enumeration, all 19 MacWilliams checks, and an independent direct enumeration of every physical error through weight five. Both are pure [[18,3,6]] codes with minimum stabilizer weight seven and D6=564. B has A7=16; C has A7=12.

The public mathematical summary is VERIFICATION_SUMMARY.md. Each derivative directory has raw rows, Pauli generators, exact transformation mask, complete enumerators, an explicit logical witness and independent verifiers. Both standalone Python verifiers need no package files or dependencies. verify_derivatives.py independently reconciles transforms, counts and text. check_clifford_phases.py gives the precise Pauli sign correction relating the all-positive Hermitian generator presentations by a global Clifford unitary.

LOCAL_UNITARY_INEQUIVALENCE.md proves that A,B,C are pairwise inequivalent under arbitrary product one-qubit unitaries and qubit permutations, using the projector weight invariant. No inequivalence under arbitrary entangling Clifford operations is claimed; the explicit transvections establish such global relationships.

All files are frozen by FROZEN_MANIFEST.json. The original A package was not changed.
