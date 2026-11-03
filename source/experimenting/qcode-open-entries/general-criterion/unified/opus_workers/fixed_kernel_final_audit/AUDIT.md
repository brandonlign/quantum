# Fixed-kernel final audit (combined prefix/suffix lex pruning)

**Verdict: SOUND for the combined-pruning argument, with the conditions and documentation changes below.** I found no correctness bugs in `common.py`, `solve.py`, `symmetry.py` or `z_enum_fixed.c`. This is a reading audit only: I reran nothing and did not reopen the `*verification.json` reports.

## Global lex-minimum argument (PROOF.md, "Combined lex conditions")

The argument is correct. Let S be the complete set of Z leaves for a fixed H, N and radical identification. Suppose:

1. every element of A and B maps S into S;
2. both tests use the same total order: lexicographic over coordinates 0..p-1, comparing integer labels.

Pick the lexicographically least tuple z* in the orbit of the group ⟨A,B⟩. Then:

- **Prefix test.** z* ≤ a·z* for every a in A. So it passes every prefix inequality, whichever site a moves, including site 0.
- **Suffix test.** If some b in B lowered z*₀, then b·z* would be smaller than z* at coordinate 0. So z*₀ is the least label in its B-orbit.
- **That least label is a kept representative.** `first_orbits` asserts that the eligible set is closed under B (`assert w in remaining`). Its representative is the orbit minimum (`root=min`).

So z* survives both prunings. You do not need a unique joint canonical form, and A and B do not need to commute. The prior audit's blanket warning is too strong when both groups preserve the same slice.

## Implementation checks against the argument

- **Witness reduction (`prefix_canonical`).** Rows before the first nonidentity row i are unchanged, so the comparison is decided at row i. Independence makes the combination ≠ xs[i]. Line 272 enforces independence by removing every span sum under `--kernel-data`. A test is skipped until its whole mask is assigned. The rejection rule `value<xs[i]` is exactly "a·xs <lex xs". The group is fully closed (`closure`), so testing all elements is correct.
- **Suffix filter.** `S^T r == r` for each radical coefficient gives r·(Sv) = r·v. So B preserves every tag test. `syndrome_orbits` asserts that Gamma and the distances are preserved. Together these preserve the eligible set and the F/tag constraints.
- **Prefix filter.** Generators are kept only if they fix each radical image modulo R. Closing only these may give a proper subgroup of the stabilizer. That is sound.
- **Identification orbits** (`identification_orbits`) are a separate, outer reduction that records every member. Applying a stabilizer subgroup inside a representative context is the standard construction.

## Required changes (documentation and proof only)

1. **PROOF.md:** retitle "(separate audit required)" and state conditions 1–2 above explicitly as the hypotheses.
2. **PROOF.md must prove, not just test, slice invariance.** Every induced action, `xs ↦ (Σ a_i[j] x_j)` for A and `x ↦ S·x` for B, must come from a physical automorphism that:
   - preserves H, R, C, N and the chosen identification;
   - maps complete (Z,Y) extensions to complete extensions.

   This is what justifies dropping all Y completions of a pruned Z leaf. It also rules out the matrix being applied transposed or inverted. A transpose is not automatically in the closed group, though the inverse is.

   Currently `verify_joint.py` only checks this empirically (`sorted(image)==range`).
3. **Independence of the evidence.** `verify_joint.py` builds its orbits from the *same* A and B used for pruning. It therefore validates the pruning logic, not whether the groups are correct or complete. Correctness of the groups rests on the asserts in `automorphisms` and `syndrome_orbits` and on `verify_symmetry.py`. Say so in the scope.

## Remaining evidence gaps (not bugs)

- Every orbit check is on small trees (84→7, 360→34, g=2 slices). No target-sized joint check has been done.
- Any production use still needs everything listed in PROOF.md's final paragraph: hash-pinned C/N classification, all H, all identification orbits, and unlimited Z and Y searches.
- `syndrome_orbits` hardcodes `d>=5` for its own representatives. `suffix_slice` does not use them, but they should not be reused elsewhere for other distances.

No production exclusion is supported or claimed by this audit.
