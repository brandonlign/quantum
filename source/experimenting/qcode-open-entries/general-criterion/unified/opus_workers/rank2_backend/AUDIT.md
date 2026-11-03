# Audit: rank2_kernel_backend/z_enum_cost.c (read-only review)

I edited no code and ran nothing heavy. These are the hashes I reviewed:

- frozen `fixed_kernel/z_enum_fixed.c` `3a81bf4a…a119`
- `rank2_kernel_backend/z_enum_cost.c` `464820c4…e790`
- binary `z_enum_cost` `552905e6…750d`
- `rank2_kernel/PREFIX.md` `e2249289…bd9cc`

The diff from the frozen source is exactly 8 hunks.

## Verdict: correct, given three caller obligations (end)

**Old mode is unchanged.** Every hunk either sits behind `use_zcost` or reduces to the original expression when it is 0. That covers the `t` formula, the `++m`, the deep-count exit, the `q>=1` bound, the first-rep threshold and the `max_t` line. Gate 1 against the frozen binary should therefore be byte-identical in stdout and in node/prune stats whenever the costs are `min(|a|,m-|a|)`.

The early exit is the one exception. With p independent reps their sum is a fifth, distinct nonzero label, so the old formula needs m deep labels. Cost mode skips the exit and simply prunes. Leaves still match, but node counts can differ when `deep_count<m`.

**Variable subset costs.** At step `depth`, `sums[depth][j]` is the XOR of `xs[i]` over the bits of `j` (the split construction puts `xs[depth]` at bit `depth`). So `j|(1<<depth)` is exactly the subset mask that includes the candidate. Each nonzero mask `a<2^p` has a unique highest bit, so it is tested exactly once, at that depth, against `d-zcost[a]`. This covers all `2^p-1` nonzero pure-Z classes.

The dependent last value is only XORed and never tested on its own: the `last<=x` check is already skipped under `--kernel-data`, and xdom/projection are forbidden. The full sum is tested once, as mask `2^p-1` with its own cost. So no forced singleton or ordering test survives.

**Threshold allocation.** `max_t = min(max(1, max zcost), d-1)`. Rows are only accessed after `distance-t<=0` has been ruled out, so `1 <= t <= min(zcost, d-1) <= max_t`. Costs `>=d` are skipped, which is correct because there is no constraint to apply. `zcost[0]=0` and `zcost[a]>=1` are enforced, so the row index `t-1>=0`. Costs above `MAX_M=10` are rejected. That is conservative: it fails loudly rather than giving a wrong answer, and 8-site costs are at most 8.

**Header.** In cost mode the header gives `p`, then `m=p+1`. The existing range check `m=3..10` gives p=2..9. The `2^p<=512` costs fit `zcost[1024]`. For q0, `bdim=2p<=13` because N<=8192, so p<=6. That includes p=6 (N=4096) and stays within the `[16]` arrays.

**First rep.** The threshold `d-zcost[1]` equals the depth-0 bound for mask 1, and `allowed` already enforces it. The validation only rejects reps that cannot be valid; completeness of the rep list is the caller's job.

**q0.** `g=t=0` is allowed only under `--z-costs`. In that case:

- `rank Gamma=bdim` must hold (full rank, exact inverse).
- All tags must be 0, and `F` has rank 0, so it is all zero.
- The kernel loop starts at `depth>=0` and imposes `<Gamma^-1 x_i, x_depth>=0` for every `i<depth`. That is exactly the isotropic prefix Gram condition.
- `global_radn=0`, so the avoid-radical leaf test is vacuous.
- Independence comes from the span removal.

For q>=1, the rank, tag, anchor, preimage and Gram checks are untouched, including the singular-Gamma `inverse_coeff` collision. The parity is well-defined there because the shifted `prev` lies in the image of Gamma.

**Forbidden flags.** `--xdom`, `--projection` and `--fixed-form` are rejected, and `--kernel-data` is required. Their R1 cost assumptions (`2+min(t,m-2-t)` and `2d-m`) are indeed invalid here.

## Caller obligations

1. **Prefix actions must preserve the cost data.** `--prefix-group` lex pruning is only sound for actions U whose coefficient transform satisfies `zcost[a]=zcost[U(a)]` for all masks, and that also preserve tags, F and the reps. With the old formula the frozen group was physical, so this held automatically. In cost mode the parent must use only physical actions that fix Z8 (PREFIX.md says the same).
2. **The first-rep list must be complete.** It must contain one representative of every orbit with `distance>=d-zcost[1]`, not `>=d-1`.
3. **Reordering the prefix basis is fine, but the costs must be re-indexed to match.** If the g independent radical-tag sites are moved first, the cost table must be permuted the same way (mask bit i means `x_i` in the new order).

## PREFIX.md math

- **Weight sum.** For `h` not in R, the four values `h_i+{I,Z,X,Y}` on each site are all different, and exactly three are non-identity. So `sum wt(h+r)=24`, and four terms that are each >=6 must all equal 6.
- **Parity.** `wt(a+b) ≡ wt a + wt b + <a,b> mod 2`, checked per site. Since every word is even, H is isotropic. This is correct.
- **q=1.** `wt(h+P8)=8-#P(h)`, and every translate has weight 6, so each letter appears exactly twice. Any two such h are related by a site permutation that fixes `<Z8,X8>`, so there is one class. This is correct.
- **Inherited claims.** The bound q<=k=3, the identity `rank Gamma=a-2q+t` and the H4 three-class / H5 one-class census come from the recorded checks. I did not re-run them.
- **Gamma dimension.** `bdim=12-q` with g=q, t=0 matches `2p-g-t` for p=6.
