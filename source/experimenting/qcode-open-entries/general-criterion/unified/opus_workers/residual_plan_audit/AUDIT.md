# Residual plan audit — verdict: SOUND

This was a read-only review. Nothing was rerun. The file hashes in `kernel_y_parity/verification.json` match the files on disk (solve, PROOF, verify, fixed_kernel/solve, common, model, two_site).

## kernel_y_parity

**Derivation.** In normal form x_i = e_i, and the DFS imposes the pair equation U_j[i] = G[i,j] ^ U_i[j] (`want=((G[i]>>j)^(v>>j))&1`). The kernel equation is XOR_{j in beta} U_j = rhs (the DFS checks `total = rhs ^ XOR values == 0`). Take coordinate i and substitute the pair equation:

sum_j beta_j G[i,j] + sum_j beta_j U_i[j] = rhs_i, so beta·U_i = rhs_i + (G beta)_i.

The j=i term is trivial because G is alternating (U_i[i] = U_i[i]). Summing over i in alpha gives beta·U_alpha = alpha·(rhs + G beta). The derivation is correct.

**Code against frozen.** The only difference from `fixed_kernel/solve.py` is the per-sector candidate filter.
- Orientation: `gb` bit i = parity(G[i] & beta) = (G beta)_i. The row indexing G[i] bit j is the same one the DFS pair mask uses, and G is symmetric anyway.
- Offsets: the candidate index v is the raw U_sector, and the offset is applied afterwards (`forward[candidate ^ offset[sector]]`). This matches the DFS lookup `allowed[sector] >> XOR values`. `rhs` is the same offset-adjusted `z ^ lift_coordinate(mask, offsets)` used in DFS `equations`. So the filter tests exactly what the DFS constraints imply.
- Indexing: after filtering, `allowed[sector]` is built from `candidate[ok]`. Before, it came from `np.flatnonzero(ok)` over `u`. This is correct because `ok` is aligned with `candidate`.
- Only the low p×p Gram block is used (beta < Q), and high rhs is still rejected by `kernel_projection`.
- An early `empty_sector` return caused by the filter is still a valid COMPLETE result, because every full solution passes the filter.

The verification results back this up: identical solution sets on all three complete map sets (one cut 920→280 nodes), and both 17,3,3 positives were FOUND with dual minimum ≥3. **No defects.**

## kernel_case_plan

**Distinguished pairs.** The coloring is `{i_p}`, other words, 24 letter vertices, 8 site vertices. This is correct for n=8: word vertices touch only letter vertices, and each letter triple touches one site vertex. So any colored isomorphism is a site permutation plus a per-site letter permutation, maps p to p', and induces a linear symplectic map of span(H). Equal (q,t,cert) therefore means the pairs (H,p) are physically equivalent. Because the map acts only on the 8 prefix sites, it extends by the identity on the suffix and preserves C and N. The q2/t0 profile indices are 4,5,6,7,9,10,11,13,14,15,16 (eleven in total), grouped into six classes: [4],[5],[6,7],[9,10,11],[13,14,15],[16]. This matches PROOF.md. All saved p words are kept (`kernel_capacity/prepare.py` iterates `row['full']`), so no p is dropped.

**Hull-3 lines.** All seven nonzero radical words are listed, and `assert len(rad)==h` guards this. The automorphisms are physical automorphisms of C (N=()), and they preserve rad. The `assert sorted(image)==range(7)` check passes. Union-find over generators gives true group orbits, and every orbit's members are recorded. Check: 5995×7 = 41965 physical lines in 41006 orbits. Assigning each orbit to its smallest-index root is fine. Nothing is excluded: only hull 5 is skipped, and that is routed to the separate q3 branches.

**Residual dependencies (inherited, not defects):**
1. The (q,t,hull,dim N) arithmetic split: hull1→{q1/t0 N=0, q2/t2 N=rad}, hull3→q2/t0 only, and dim N = 1 for hull3/q2. This is asserted, not re-derived here.
2. Completeness of `kernel_capacity/profiles.json`, which relies on the parents `gauge_prefix_candidates.json`, `cells_d5_7_3.json`, and `all_critical_dim4.json`, and on the q1 classes already being R-preserving representatives.
3. The suffix9 hull-1/3 class census (26752 / 5995).
4. Downstream extension search must treat H and N as subspaces, so its result does not depend on the basis chosen.

**Minor (non-blocking):**
- `plan.json` hashes `profiles.json` but does not check it against the hashes `profiles.json` itself records, and it omits the upstream parent files and `kernel_case_plan/PROOF.md`. Consider adding these for provenance.
- PROOF.md should state dependency 1 explicitly.

As both documents say, neither artifact proves any target exclusion. Each representative still needs a complete exact extension search over all its N orbits.
