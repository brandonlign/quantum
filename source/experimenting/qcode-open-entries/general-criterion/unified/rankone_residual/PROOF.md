# Residual q1/q2 exact fixed-H,N search: coverage argument

Scope: exclusive rank-one column R=<Z^8>, cases q1/t0/hull1 (N=0),
q2/t2/hull1 (N = the unique radical line) and q2/t0/hull3 (N = each of the
seven radical lines, reduced to physical Aut(C) orbits). This namespace does
not by itself exclude [[17,3,6]]. LIMIT/ERROR/missing never exclude.

## Inputs (all hash-pinned in every record)

- `kernel_case_plan/plan.json`: every suffix class of hull 1/3 with its
  N-orbits (members recorded) and the distinguished (H,p) classes:
  4 q1/t0, 2 q2/t2, 6 q2/t0. Equal distinguished certificates are physically
  equivalent pairs (kernel_case_plan/PROOF.md; residual plan audit SOUND).
- `kernel_case_plan/provenance.json` (PASS, 25 nested digests) validated at
  start of every run and aggregation.
- Z backend `kernel_z_compact/z_enum_compact` (stream/count identical to the
  frozen fixed_kernel backend; verification PASS). Y solver
  `kernel_y_cache/solve.py` = `kernel_y_parity/solve.py` with an LRU on the
  immutable standard(m) frame (verification PASS).
- Prefix contexts / identification orbits: `kernel_context_cache/common.py`
  caches the prefix-only calculations from `fixed_kernel/common.py` and
  substitutes the actual suffix radical basis in its original order. Its
  proof and 24 exact returned-data/orbit comparisons are recorded separately.
  Suffix first-label slices still use `fixed_kernel/symmetry.py` unchanged.

## Task decomposition and why it covers every code

Two native backends are registered and hash-pinned: compact and compact with
the exact four-free MITM feasibility test. Every record identifies its backend.
Both enumerate the same complete Z set. The controls compare complete LEAF
streams from both on all pinned positive families; the separate joint checks
also cover artificial distance-six singleton/pair thresholds. Complete chunks
from either registered backend may be reused together, since their coverage
and chunk partition are identical. `joint-auto` selects MITM only for profiles
0 and 8, where paired timing probes showed a clear gain. Other profiles retain
compact-only. Changing this runtime selection does not change the semantic
source fingerprint or turn an unfinished chunk into an exclusion.

For each suffix class C, branch (q,t), N-orbit representative, distinguished
H, we compute `kernel_contexts` = ALL GL(g,2) radical identifications
(1, 1, 6 for g=0,0,2... specifically g=q-t: q1/t0 g=1 -> 1, q2/t2 g=0 -> 1,
q2/t0 g=2 -> 6). `identification_orbits` merges them under the full physical
Aut(H,p) (all members recorded; the aggregator asserts the union equals all
contexts). A complete extension in identification j is transported by a
physical prefix map to one in the orbit representative, preserving C, N
(the map is identity on the suffix), distances and exact kernel
(fixed_kernel/PROOF.md, "Prefix actions"). So searching representatives only
is complete.

Within a context, prefix lex pruning (fixing each identified class mod R) and
minimum first-label representatives of the radical-fixing Aut(C,N) subgroup
are combined by the joint lex-least argument in fixed_kernel/PROOF.md
("Combined lex conditions"), exactly as in the production hull3 runner.

The first-label representative list is split into consecutive chunks of
CHUNK=64. The Z tree restricted to a subset of first labels is exactly the
union of the restricted trees (the first-orbit mask acts only at depth 0),
so the slice is exhausted iff every chunk is COMPLETE with every Y leaf
COMPLETE and `fixed_z_excluded`. The aggregator recomputes the chunk count
and the digest of the eligible set / orbits from scratch.

Exactness of K and N, Y completion, M(H)=0 and all distance conditions are
those of fixed_kernel/PROOF.md and kernel_y_parity/PROOF.md.

## Witness path

Any FOUND first writes `candidates/CANDIDATE_<task>.json` atomically with
the raw columns, before any check. Then `physical_witness` independently:
constructs S,D, requires rank S=14, rank D=20, rank M=14-q (exact N), S
isotropic, S commuting with D, and enumerates all 2^20 words of D to get the
minimum weight outside S. A physically valid distance-six code is a target
witness regardless of profiles. Errors in this check are recorded and the
candidate file is retained; the runner stops on any FOUND.

## Controls (`controls.py` -> `controls.json`)

For two hull-1 and two hull-3 classes, every branch and every distinguished
H: build a genuine [[17,3,3]] code with that exact C, N, H by an explicit
presymplectic isometry; locate its identification; require the production
compact Z backend (no group, first 7/6/5 labels pinned) to emit its Z tuple;
require the production Y solver to find a completion; physically verify both
with the 2^20 dual check at distance three; and replay all identification
and first-label orbit coverage. These are distance-three controls, not target
constructions. Production `run.py` refuses to start unless controls PASS.
