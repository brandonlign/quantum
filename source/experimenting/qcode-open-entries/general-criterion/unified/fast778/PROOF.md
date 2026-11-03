# S6 coverage for the 778 frame

The frame R=<Z^8,X^6 Z I> is fixed by every permutation of its first six
sites. Its adapted pure-Z basis is exactly the same as in the 688 frame:
Z6, Z0+Z5,...,Z4+Z5. The increasing-differences and partial-anchor argument
in speedup_review/PROOF.md therefore applies unchanged, with the same fixed
first image and the same disabled necessary-only unary Y filter.

The existing Stage-2 minimum outside R is eight. Since R has exactly two
weight-seven words and the global bound A6+A7<=2 holds, this also fixes the
global profile (0,2). Requiring exact final outside counts (0,0) does not
strengthen the existing mu_out=8 condition. No changed purity assumption
or inherited minimum-eight full-plane exclusion is introduced.

All 52 retained classes and all 650 original first roots are preserved in
the new manifest. Their original physical inputs and the original frozen
source hashes are copied/read without altering Stage 2. All new results
have a separate namespace. Exhaustion and Y-completion evidence use the
same verified streaming/checkpoint pipeline as the isolated 688 runner.

Physical S6 generators are verified on every quotient coset and all four
R-translate weights. A physical n=14 graph survives nontrivial transport
and exact profile counts. The small fixture worker reconstructs genuine
solutions, and rejects forged exclusion records. No completed target
exclusion is supplied by this implementation validation alone.
