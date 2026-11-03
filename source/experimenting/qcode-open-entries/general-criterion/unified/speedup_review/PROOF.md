# Coverage of the isolated S6-reduced 688 search

For R=<Z^8,X^6 I^2>, permutations of sites 0,...,5 fix both generators.
They preserve the physical symplectic form, weight, R, and the full pure-Z
quotient Lagrangian. The adapted pure-Z basis is
Z6, Z0+Z5, ..., Z4+Z5. Its first image X0 is fixed. The other five images
are the differences from an anchor in a set of six affinely independent
points. All six Z images are independent in a hypothetical extension:
every nonzero pure-Z quotient class has cost at most four, below distance six,
so its syndrome cannot vanish. The original span-exclusion rule is retained.

Choose the lexicographically least tuple of five ordered differences over
this physical S6 orbit. It is increasing. At each depth, the new backend
enforces this increasing order for X1,...,X5. For any assigned anchor,
it sorts the assigned other points after translation and rejects if that
tuple is smaller than the assigned current tuple. Unknown future points
can only lower the sorted translated prefix, because adding elements to
a sorted list cannot increase its first entries. Thus a witnessed smaller
assigned prefix excludes every full extension of that prefix from being
the chosen orbit minimum. At full depth this is the complete anchor test.

The first image remains fixed throughout. Therefore every existing eligible
first-image root retains a representative of every physical solution in
that root; the roots are still disjoint and exhaustive. No further suffix
symmetry is introduced. The exact outside weight profile is invariant under
physical site permutations. The full Y images are transported with the
entire symplectic quotient action, including any Z/Y shear. The exact frozen
Y completion keeps that freedom and every complementary sector.

Selected unary Y domains of the predecessor are not invariant as relaxed
Z constraints under arbitrary basis changes. The new canonical lane disables
this optional necessary-only filter. The predecessor's other necessary tests
remain necessary for every full physical extension in its transformed basis.
This proves preservation of full solutions. Equality of raw relaxed Z-leaf
sets with the predecessor is not required or asserted. Additional necessary
leaves must still be completed exhaustively by the frozen Y solver.

The copied Stage-3 worker requires a COMPLETE enumerator terminal summary,
unique matching leaf/completion hashes, complete unlimited Y completion,
physical reconstruction of every returned solution, and exact final profile
equality. A root is excluded only if this pipeline exhausts with zero accepted
solutions. Aggregation requires every one of the 1,950 unchanged first roots.
LIMIT, ERROR, missing summaries and unfinished roots remain unresolved.

All new sources, binary and copied inputs live outside mixed_frames. Both
the predecessor implementation fingerprints and new implementation hashes
are frozen in the new manifest. Its namespace does not reuse predecessor
checkpoints; current original searches retain their original provenance.

Validation includes complete small enumeration and constructive physical
positive controls, disabled-symmetry equality against the frozen binary,
all 20,480 generator/quotient checks, and the actual worker/ack/cache pipeline.
These are self-checks and a mathematical argument, not an independent audit
or a new full [[17,3,6]] exclusion.
