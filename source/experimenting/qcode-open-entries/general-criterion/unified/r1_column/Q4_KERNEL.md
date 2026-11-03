# Kernel dimension four is impossible in the rank-one weight-eight frame

13 October 2026. Found by a Sol 6.1 worker (sol_workers/C_rank_one_attack.md, §3) and re-derived and checked by hand. It is a small structural lemma, not an exclusion of [[17,3,6]].

**Claim.** Let S be a stabilizer code with parameters [[17,3,≥6]], let T be an eight-site support, and suppose S ∩ V_T = R = ⟨Z^8⟩ *exactly*. Then the extension kernel has dimension q ≤ 3. KERNEL_BOUND.md's Griesmer bound gives only q ≤ 4.

Purity is not assumed. Neither is any bound on the minimum stabilizer weight, so the lemma covers the minimum-six weight-eight frame of WEIGHT8_RECENTER.md.

**Proof.** Let H = S^⊥ ∩ V_T. Then H ∩ S = R, H/R ≅ ker M, and so dim H = q+1. Every word in H \ R is a nontrivial logical operator, so it has weight at least 6. Suppose q = 4, so |H| = 32.

1. *Averaging.* In an additive code over GF(4), each coordinate that is not identically zero is nonzero in exactly 3/4 of the words. So Σ_{h∈H} wt(h) ≤ 8·32·3/4 = 192.
2. *Lower bound.* The 30 words of H outside R have weight at least 6. Each odd-weight word has weight at least 7, and Z^8 has weight 8. If O is the number of odd-weight words, Σ wt ≥ 30·6 + 8 + O = 188 + O, so O ≤ 4.
3. *Parity.* The map h ↦ wt(h) mod 2 is a binary quadratic form whose polar form is the symplectic product. This is checked for all 16 pairs of single-site Paulis. A nonzero quadratic form on F₂⁵ takes the value 1 at least 8 times; an exhaustive check of all 2^15 forms confirms this. So the form vanishes and H is isotropic.
4. Because H ∩ S = R, H/R embeds as a 4-dimensional isotropic subspace of the 6-dimensional logical symplectic space N(S)/S. But a maximal isotropic subspace there has dimension 3 = k, a contradiction. ∎

The finite ingredients are checked by `verify_q4.py`, which writes `verify_q4.json` with status PASS.

**Effect.** The minimum-six arithmetic kernel screen admitted (q,t,h,dim N) = (4,2,5,3) and (4,4,3,3). Both are now removed. In the recentred rank-one hull-five case this forces q = 3, t = 0 and dim N = 2, so the 2,048-label fixed-N reduction applies again. Its N words only need weight ≥ 6, though, so the old heavy-weight-eight plane coverage does **not** transfer. `kernel_bound.py` has not been changed.
