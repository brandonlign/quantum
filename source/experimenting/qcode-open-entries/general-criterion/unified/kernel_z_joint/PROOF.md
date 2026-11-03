# Compact Z columns and full-label four-free feasibility

This experimental composition retains the compact backend's singleton domain
H={v:delta(v)>=d-1} and the MITM backend's necessary depth-four feasibility
test, with all original recursion tests still present. It is not production
and proves no exclusion without a validated complete search.

Every free singleton label lies in H, so the four singleton domain bitsets
may use compact columns. Every bit is converted to its actual syndrome label
before XOR, Gamma, tag, pair or bucket operations. Pinned prefixes and the
brute-force implied fourth label use the inverse label-to-column map. Padding
bits are masked before domain filtering.

Pair sums are different: their threshold is d-2, so they need not lie in H.
Their predicates operate on FULL syndrome labels. `pair_allowed(v)` directly
checks every assigned subset A: delta(v+sum(A)) >= d-min(|A|+2,8-|A|-2).
These are exactly the original MITM pair-mask predicates. Epoch-stamped arrays
cache those decisions within one node and mark membership in the full-label
assigned span. Caches never carry predicates to a different node, and epoch
wrap clears the stamps. Skips and work overflows still return feasible.

Consequently each original MITM domain, span test and pair test is identical
as a set of full labels. All matching, tags and bilinear pair tests are copied
unchanged from the validated MITM backend. Complete leaf streams must match
both it and the unoptimized backend; node counts match the original MITM
backend when the same cap and node limit are used. Timing needs separate
paired measurements because candidate-domain compression and lazy pair tests
have different overhead.
