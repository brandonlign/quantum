# Seventeen-qubit computation: map and rank-one checker

## 1. The computation behind the [[17,3,6]] nonexistence proof

Original evidence lives in
`experimenting/qcode-open-entries/general-criterion/unified/` (paths below are
relative to it). Each row is checked by `verify.py` at the workspace root.

| Step | Input / record | Recorded result |
|---|---|---|
| Purity and light words (Prop. 5.1) | `verification/certificates` | 7 certificates pass; A6 + A7 ≤ 2 |
| Shortened-code classification | `parent_mass_check/chain.json`; `opus_parallel/suffix9_dimension5/suffix9_dim5.json`; `suffix9_classes.json` | 32,865 dimension-5 classes (hull 1: 26,752; hull 3: 5,995; hull 5: 118); 56 dimension-6 classes |
| Marked local kernel types | `prefix_mass_check/verification.json` | PASS |
| Rank two, profile (6,6,8) | `mixed_frames/stage2/` | 48 classes, 641 roots, 0 Z leaves |
| Rank two, profile (7,7,8) | `fast778/` | 52 cases, 650 roots |
| Rank two, profiles (6,8,8) | `speedup_review/` | 156 cases, 1,950 roots |
| Rank two, full plane (8,8,8) | `rank2_kernel/combined_888.json` | 56 classes, 1,375 slices, 2,960 Y leaves rejected, 72,827 hull-0 roots |
| Rank one, (q,t,h) = (3,0,5) | `radical_projection/aggregate.json` | 118 classes; 5,749 screened, 7,065 searched |
| Rank one, (3,2,3) | `fixed_kernel/hull3_aggregate.json` | 5,995 classes; 20,791 screened, 20,215 searched |
| Rank one, residual cases | `oracle_residual/aggregate.json` (records in `evidence/search_records.sqlite`) | 32,747 classes, 2,419,796 records |

| Replay of the saved evidence | Result | Time |
|---|---|---|
| `verify.py` (all branches re-aggregated; 8 residual classes, 1,059 records, replayed) | PASS (`provenance/portable_replay.json`) | 71 s |
| `verify_all_record_hashes.py` | PASS (`provenance/all_record_integrity.json`) | 171 s |
| `verify.py --deep` (every residual record) | started 16 Nov 2026; log `provenance/deep_replay_20261116/deep_replay.log`; result goes to `provenance/portable_deep_replay.json` | running |

## 2. Exact rank-one test without kernel cases

`rankone_exact.py` (standard library + numpy) decides, for one shortened code C,
whether a rank-one extension with R = ⟨Z^8⟩ exists. It uses only Theorem 2.1
(extension criterion), Proposition 5.1 (purity) and the classification of C.
It needs no kernel cases, local types, recentering or saved records.

Basis of E = R^⊥_T / R: z_i = Z_i, x_i = X_i X_7 (i = 0..6). Every class is
X_B Z_A with B even, and cost(X_B Z_A) = |B| + min(|A∖B|, |B^c∖A|).
δ is computed by breadth-first search over the 8,192 syndromes.

* **Stage 1 (Z-half).** ξ_i = M(Z_i), with Σ ξ_i = 0 and
  δ(Σ_{i∈A} ξ_i) ≥ 6 − min(|A|, 8−|A|) for every proper nonempty A
  (thresholds 5, 4, 3, 2 for |A| = 1..4). The images are sorted, using the
  symmetry under site permutations of T. Γ = MJM^T forces the annihilator U of
  span(ξ) to be Γ-isotropic. After k images this gives
  rank(Γ|U_k) ≤ 2(7−k), and when equality holds the later images lie in the
  annihilator of rad(Γ|U_k).
* **Stage 2 (X-half).** η_a = M(x_a) = base_a + Ξ(Y_a) with
  (Y_a)_b + (Y_b)_a = f_a^T Γ f_b. The |B| = 2, 4 distance conditions become
  domain conditions on partial sums of the Y_a, and an exhaustive backtrack
  over Y_0..Y_6 solves them. Conditions with |B| ≥ 6 are automatic.

Together the two stages are equivalent to the Theorem 2.1 condition.
`full_check` re-verifies any returned M directly: the Gram equation and all
2^14 distance conditions.

## 3. Results

| Program | What it did | Result | Time |
|---|---|---|---|
| `survey_rankone.py` | For all 32,865 classes: deep holes (δ ≥ 5), the graph δ(u+v) ≥ 4, and an 8-tuple meeting every Stage-1 subset condition, ignoring Γ | deep holes 476–1,702 (mean 858.5); edge density 0.497–0.631; a tuple was found for all 32,865 classes (class 24383, hull 5, in 10.3 s without a time cap: `circuit_24383.log`) | 905 s, 5 workers |
| `estimate_stage1.py` | Random-probe estimate of the full Stage-1 tree (d = 6, no kernel restriction), 200 probes on 60 + 30 + 10 classes | mean nodes per class: hull 1 3.3·10^8, hull 3 4.2·10^8, hull 5 8.5·10^9 (max 4.9·10^10) | — |
| `control_d5.py` | Built a pure [[17,3,5]] code from a d = 5 extension (brute force: rank 14, isotropic, all 2^20 logical-space words weight ≥ 5), scrambled it, and for every rank-one weight-8 stabilizer word recovered the true M | 131 words (hull 1, d(C)=5: 107; hull 3, d(C)=5: 23; hull 1, d(C)=6: 1); true ξ and η recovered by Stage 1 and Stage 2 in every case: PASS | — |
| `control_dfs.py` | Stage-1 completions of a fixed 5-prefix (d = 5) vs. brute force over all (ξ_5, ξ_6) with every subset condition and exact isotropy | classes 0, 1, 5, 150: 1,549 / 2,215 / 1,932 / 705 completions, identical sets: PASS | ~130 s |
| `run_exact.py 150 --ndim 2 --nindex 0,1` | Exact test of class 150 (hull 5) with im M = N^0, for 2 of its 155 planes N ⊆ rad(P) | 0 Z-halves, 0 extensions for both planes; ~1.25·10^6 Stage-1 nodes each | 165–269 s each |

**What transfers from the 16-qubit proof.** The pair/circuit view of the
Z-images transfers as a correct necessary condition, but at m = 8, d = 6 it
excludes none of the 32,865 classes: every class admits an 8-tuple meeting all
Z-only subset conditions. What makes the problem decidable is the Gram
(isotropy) condition plus Stage 2. Together they give an exact per-class test
that replaces the kernel-case split.

**Cost.** Summed over all classes, the estimated Stage-1 tree has about
1.2·10^13 nodes. At 4,600–7,600 nodes/s in Python this cannot be finished here.
Restricting to a kernel case shrinks the tree a great deal. Class 150 used
about 1.25·10^6 nodes per plane. Without the restriction, the ten sampled
hull-5 classes have estimated trees of 3.0·10^8–4.9·10^10 nodes.

## 4. Remaining gaps for an independent rank-one check

* The exact test has been completed for class 150 on 2 of its 155 planes,
  not for any class in full and not for the other 32,864 classes.
* An unrestricted full run needs about 1.2·10^13 Stage-1 nodes; a run over
  kernel cases needs the N-orbit representatives of the paper's split as
  input (`--ndim`, `--nindex`).
* The light-word bound A6 + A7 ≤ 2 (Prop. 5.1) is not used as a filter.
* Rank-two extensions (profiles (6,6,8), (7,7,8), (6,8,8), (8,8,8)) are
  outside this checker.

## Files

`rankone_exact.py` (core), `run_exact.py` → `run_exact.jsonl`, `run_exact.log`;
`survey_rankone.py` → `survey_rankone.json`, `.log`; `circuit_24383.log`;
`estimate_stage1.py` → `estimate_stage1.json`, `.log`;
`control_d5.py` → `control_d5.json`, `.log`; `control_dfs.py` → `control_dfs.json`, `.log`.
