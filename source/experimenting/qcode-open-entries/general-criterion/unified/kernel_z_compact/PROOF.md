# Compact candidate columns, unchanged syndrome algebra

This backend is a representation-only copy of fixed_kernel/z_enum_fixed.c.
It applies ONLY to that backend's physical cost-one singleton Z generators,
not the variable-cost rank2 backend. Every chosen X_i must satisfy
delta(X_i)>=d-1, by the subset test for its singleton. Thus it is exact to
restrict candidate columns to the sorted list H={s:delta(s)>=d-1}.

Every distance and parity bitset now has one bit per H member. Its ROW index
remains a full syndrome label. XOR sums, Gamma coefficients, radical tags,
prefix action coefficients, chosen labels and output labels all retain their
original values. The candidate-bit index is converted to its label before
any algebra. Span exclusions, first masks and pinned prefixes use the inverse
label-to-column lookup. Ordered enumeration uses upper_bound on actual labels.
Unused padding bits in the final word are masked at every node.

For each row (t,s), candidate bit j is exactly
delta(H[j] XOR s)>=d-t. The Gamma parity bit is exactly c dot H[j]=0.
These are the original predicates restricted to a subset that every original
candidate already satisfied. Sorted H preserves traversal order. Consequently
the LEAF stream, visited recursion nodes and all combinatorial prune counters
are identical to the original for the same input and NODE limits.

Only the old xdom option is rejected: its complementary-Y domains range over
all labels and cannot share the compressed Z rows. It was measured slower and
is not used by the fixed-kernel runner. Projection tags, when requested, still
use actual full syndrome labels and therefore remain unchanged.

The representation reduces bitset words and precomputation from N to |H|
columns. For target R1 suffix classes, |H| is usually about 500–1600, while
N is 8192 before the exact-N quotient. No speedup is claimed until paired
complete-tree measurements pass.
