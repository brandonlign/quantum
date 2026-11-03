# Audit: suffix_mass_check as a replacement for a second direct d6 census

**Verdict: SOUND** relative to the residual dependencies listed in §8. A complete labelled-mass check of the (8, dim5, d5) parents, followed by full new-certificate coverage and per-class multiplicity equality for every valid canonical lift, is enough to prove the 32865-class (9, dim5, d6) cell is complete and has no duplicates. A second direct d6 census is not logically needed. Nothing below excludes any quantum target or makes a novelty claim.

What I read: PROOF.md, common.py, verify8.py, verify9.py, verification8.json, verification9.json (summary fields), and release16/census/check_cells.{py,log} (only to pin down the dependency). I re-ran nothing heavy.

Recorded outputs:
- verification8.json: PASS_RELATIVE_TO_CERTIFIED_RELEASE_PARENTS. LHS = RHS = 1593099450132480. 26299 distinct children. The dim3 check found 116280 multisets, 2171 feasible, 32 classes, and the certificate sets are equal. Run time 2.78 s.
- verification9.json: PASS. 26299 parents, 1791031 ordered pairs, 299175 canonical lifts, 32865 source classes. Multiplicities sum to 299175, the minimum is 1 (so every class is hit), and the maximum is 88. The hull counts are 825 + 53957 + 244393 = 299175. Run time 34.7 s.

## 1. Shortened-child recurrence (verify8)

Take a child C of type (8,5,≥5). Let S be the subcode with letter I at site 8, restricted to sites 1–7, and let V = π₈(C) ⊆ F₂², with r = dim V ∈ {0,1,2}. Then:
- dim S = 5 − r, and S has distance ≥ 5.
- The words of C with last letter v form a coset ψ(v) of S. The map ψ: V → F₂¹⁴/S is linear and injective.
- C has distance ≥ 5 exactly when every nonzero ψ(v) has minimum weight ≥ 4 (the "deep" condition).

The map C ↔ (S, V, ψ) is a bijection on labelled subspaces. Each C has exactly one shortening at site 8. The number of extensions does not change under equivalence of S, so it is enough to compute it once per orbit representative. ✔

The quotient is implemented correctly. Labels are symplectic syndromes against a basis of S^⊥s, which is a surjective map F₂¹⁴ → F₂^{9+r} with kernel S. The code asserts this kernel equals span(g). All 4⁷ words are assigned to a label, so every coset minimum is exact. ✔

## 2. The factors 3 and 6, and the Fourier cubic count

- r=0: there is one extension. ✔
- r=1: there are 3 choices of the line V and #deep choices of ψ(V), giving 3·#deep. ✔
- r=2: ψ is determined by a=ψ(X) and b=ψ(Z), with ψ(Y)=a+b. The count is #{(a,b): a, b, a+b all deep}. If a=b then a+b=0, which is never deep, so injectivity is automatic. This equals 6 × (number of deep planes). The code asserts count%6==0. ✔
- The Fourier identity Σ_{a,b} D(a)D(b)D(a+b) = |Q|⁻¹ Σ_u D̂(u)³ uses the unnormalised WHT, and fwht is a standard butterfly. The code checks it exactly against a direct ordered-pair count and checks divisibility. ✔ Because the direct count is also asserted, an error in the Fourier step cannot pass silently.
- The orbit weight 7!·6⁷/|Aut(P)| is correct for G₇ = S₇ ⋉ S₃⁷ acting faithfully on subspaces of F₂¹⁴. ✔

So LHS is the exact number of labelled (8,5,≥5) subspaces, provided the parent list is complete with one entry per class (§4).

**Completeness logic.** RHS is a sum of orbit sizes over children that are valid and pairwise inequivalent (certificate not already seen), so RHS ≤ the true total. If LHS = the true total, then LHS = RHS forces the child list to contain every orbit. ✔

## 3. Star-graph encoding, distinctness and Aut

The graph has codeword vertices, 3n letter vertices and n site vertices, each kind with its own colour.
- Each letter vertex is adjacent to exactly one site vertex. So a colour-preserving automorphism gives a site permutation σ together with a bijection from site j's letters to site σ(j)'s letters. That is an element of S_n ≀ S₃, and S₃ = GL(2,2) = Sp(2,2), so these are exactly the physical local operations.
- Distinct codewords have distinct letter neighbourhoods. The zero word is the unique isolated codeword vertex. So the map Aut(graph) → G is injective, and its image is exactly Stab_G(C). Idle sites and unused letters are handled correctly: the S₃ on an idle site still appears.
- Equal certificates therefore mean physically equivalent codes, and the graph group order equals the physical |Aut|. ✔
- This encoding is independent of the old triangle encoding. The orders are also checked against the saved orders and divide n!·6ⁿ.

Minor: `round(grpsize1*10**grpsize2)` uses floating point. Group orders here are at most about 7·10¹⁰, so exactness is very likely. Any error would also break the exact integer mass equality or the divisibility assert, so this cannot cause a false PASS.

## 4. Parent completeness (lengths 7, dims 5/4/3)

- **dims 4 and 5:** Their certificate sets are asserted equal to those of release16 cells.json '7,4' and '7,5'.
  - check_cells.py confirms that '7,4' has 43 distinct representatives equal, as a set, to its mass-checked census. Since the list here also has 43 entries, it has no duplicates.
  - '7,5' is not included in the set comparison inside check_cells.py. However, the census log shows n=7 k=5 has 1 class with the mass identity holding. verify8 asserts that the single parent is valid (rank 5, d ≥ 5), so it must be that class.
  - **The 7,4 and 7,5 cells are the only inherited completeness.**
- **dim 3, independent:** Each column of a 3×7 generator is a pair of forms on F₂³. Up to local GL(2,2) it depends only on the row space, giving 1 + 7 + 7 = 15 types. Each plane is represented uniquely by a<b<a^b.
  - Any code with any basis gives a multiset of types, and the code built from those types is physically equivalent to it. So the C(21,7) = 116280 multisets are exhaustive.
  - The weight prefilter is exact: weight(c) is the number of columns that hit c. A rank defect would create a weight-0 nonzero message, so it is filtered out automatically.
  - The 32 new certificates equal the saved set, and since the list has 32 entries, it has no duplicates. ✔

## 5. Puncture injectivity and tail-form coverage (verify9)

- **Puncture:** If C has type (9,5,≥6), then a nonzero kernel of puncturing site 9 would have weight ≤ 1. So the punctured code P has dimension 5 and distance ≥ 5, which makes it equivalent to some stored parent.
- **Lift:** Apply that equivalence to sites 1–8. Then C = {(p, f(p))} for some linear f on the stored basis, which is a pair of forms (a,b) on F₂⁵. All 1024 pairs are counted. Up to GL(2,2) at site 9 there are 1 + 31 + 155 = 187 row-space representatives, with the same unique-representative rule as above. ✔
- **Validity:** The lift has distance ≥ 6 exactly when f(p) ≠ 0 for every p of weight 5. Since the parent has distance ≥ 5, smaller weights cannot occur. The bad/zero masks are indexed consistently with xor_table (bit i of c corresponds to g_i), and so is the bit layout of the lift (x at bit 8, z at bit 17). Distance ≥ 6 is re-asserted with span. ✔
- **Conclusion:** Every class occurs among the lifts. Each lift's new certificate must be in the source set, and every source class must be hit, with all source representatives checked as valid and distinct. So the cell has exactly 32865 classes. The per-class multiplicity equality and the hull checks are extra consistency evidence; completeness does not depend on them.

## 6. Circularity hazards

| Hazard | Status |
|---|---|
| verify8 uses the old child list | Only as a target. Validity, distinctness, Aut and mass are all recomputed, and the mass argument cannot be satisfied by an incomplete list. Not circular. |
| verify9 uses the old 32865-class list | Only as a target. Both inclusions are asserted with new certificates. Not circular. |
| Parent files for dims 4/5 come from opus_parallel | They are tied to release16 by certificate-set equality. Their own `mass_equal` flags are not relied on. |
| Saved labelled_ext and aut_orders | These are extra equality checks only. |
| verify9 using a different 8-site list than verify8 checked | Prevented. verify9 re-hashes every file listed in verification8.json (including cell_8_5.json, common.py and model.py) and freezes the source file's hash before and after the run. |
| Shared code with the old pipeline | model.py (span, reduced_basis, dual, gram, weight), xor_table and pynauty are common-mode. Partial cross-checks: the syndrome kernel equals span; xor_table and span are both used in verify9; the graph is new. Low residual risk. |
| Same method family as release16 | release16 is itself a lengthening/mass census. Its *code* is independent here, but the 7,4 and 7,5 trust is inherited, not re-derived. |

## 7. Optional hardening (not blocking)

1. In verify9, assert Σ over valid forms of (1 for zero, 3 for line, 6 for plane) = ordered_count (1791031). This ties the 187-form reduction to the 1024-pair count.
2. In verify8, assert directly that the parent certificates within each r are distinct. This is currently implied by list length and set equality.
3. To remove the release16 dependency, re-derive 7,4 and 7,5 independently. Options: the 51-type multiset method for dim 4 (C(57,7) ≈ 2.6·10⁸, not run here because heavy searches are excluded), or a verify8-style shortening mass from the 6-site cells.
4. Wording in PROOF.md: say explicitly that after transporting the equivalence to the stored generator basis, every linear tail is a pair of forms relative to *that* basis. Also note that '7,5' is accepted through validity plus the one-class mass log, not through the set comparison in check_cells.

## 8. Exact residual dependencies

1. **release16 completeness of (7, dim4, d5) [43 classes] and (7, dim5, d5) [1 class].** This covers release16/census/mass_census.py and its own cell chain, as pinned by the cells.json hash 181b2754…
2. Correct pynauty 2.8.8.1 `certificate` and `autgrp` on coloured graphs.
3. Correct model.py helpers (span, reduced_basis, dual, symp, gram, weight), mixed_frames.frames.xor_table, numpy, and Python integer arithmetic.
4. The integrity of the hash-pinned data files: cell_8_5.json, cells_d5_7_{3,4,5}.json, suffix9_dim5.json.
5. The equivalence notion is S_n ⋉ S₃ⁿ acting on additive codes (min weight over all nonzero codewords). The result says nothing about quantum targets beyond this cell.

Independent of these: the 7-site dim-3 parents, every coset minimum, the r=1 and r=2 extension counts, the 8-site mass identity, the new certificates and Aut orders, and full coverage of all valid tail forms.
