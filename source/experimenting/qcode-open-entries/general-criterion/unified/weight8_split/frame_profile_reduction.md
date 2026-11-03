# Exact profile/frame reduction for pure [[17,3,6]]

Exploratory computer-assisted result, 13 October 2026. This is not a nonexistence proof.

Using the exact rational Farkas verifier in `split/system.py`, I fixed the pure weight enumerator entries A1=...=A5=0 and enumerated the six possible pairs `(A6,A7)` with A6,A7≤2. For each pair I added one of the five local rank frames forced by a weight-eight word:

- r=1: `ZZZZZZZZ`
- r=2: `ZZZZZZZZ, XXXXZZII` (weights 6,6,8)
- r=2: `ZZZZZZZZ, XXXXXXII` (weights 6,8,8)
- r=2: `ZZZZZZZZ, XXXXXXZI` (weights 7,7,8)
- r=2: `ZZZZZZZZ, XXXXXXXX` (weights 8,8,8)

`CERT` means the split-enumerator necessary relaxation is exactly infeasible. It is a valid exclusion of that frame/profile combination, not necessarily of a quantum code, because the relaxation omits some full extension conditions.

| `(A6,A7)` | r=1 | (6,6,8) | (6,8,8) | (7,7,8) | (8,8,8) |
|---|---:|---:|---:|---:|---:|
| (0,0) | open | CERT | CERT | CERT | open |
| (0,1) | open | CERT | CERT | CERT | open |
| (0,2) | open | CERT | CERT | open | open |
| (1,0) | open | CERT | open | CERT | open |
| (1,1) | open | CERT | open | CERT | open |
| (2,0) | open | open | open | CERT | open |

Thus the current exact relaxation reduces the possible rank-two frame choices as follows:

- `(0,0)` and `(0,1)`: only (8,8,8) remains among the rank-two choices;
- `(0,2)`: (7,7,8) or (8,8,8);
- `(1,0)` and `(1,1)`: (6,8,8) or (8,8,8);
- `(2,0)`: any of the (6,6,8), (6,8,8), or (8,8,8) frames remains; (7,7,8) is excluded.

A separate exact run proves the global pure enumerator relaxation infeasible for `A6≥3`, `A7≥3`, and `A8≤119`. Therefore the six pairs above are exhaustive at the ordinary enumerator level, subject to the exact certificate scope.

## What this does not establish

The r=1 column is not closed. The rank-two results do not apply to the whole minimum-6, minimum-7, or minimum-8 branch; they apply only after selecting a weight-eight word and a specified local rank frame. The output records Farkas status but not a separate certificate-file replay yet. Before citing this as final proof evidence, dump the exact multipliers and run an independent reconstruction checker.

## Next exact work

1. Dump multipliers for the 11 `CERT` cells and independently rebuild each split system from labels.
2. Extend the cost-table/extension search to the mixed frames `(6,6,8)`, `(6,8,8)`, `(7,7,8)`.
3. For (8,8,8), use the profile reduction to prioritize `(0,0)`, `(0,1)`, then the remaining `(0,2),(1,0),(1,1),(2,0)` profiles.

## Correction and independent replay (13 October, Sol 6.1 audit + root check)

Source: `../sol_workers/A_certificate_audit.md`.

- **The table has 13 CERT cells, not 11.** Each one is a plain local word-count obstruction. The frame forces A6 ≥ 2 for (6,6,8), A6 ≥ 1 for (6,8,8), and A7 ≥ 2 for (7,7,8). A two-row certificate, unit·(−h) + weight_w·1, proves each cell; it needs no MacWilliams or shadow rows. All 13 pass the unchanged standalone checker `pure/check_general.py`. The bundle is in `frame_count_certificates.json` and the PASS lines are in `frame_count_certificates.log`.
- **Exhaustiveness of the six profiles.** The bounds A6 ≤ 2 and A7 ≤ 2 alone allow nine pairs. The exact identity A6 + A7 ≤ 2 closes the gap. Its multipliers are unit −1835008, total 116, low rows (60, 25, 11, 3, 1), weight rows (−851968, −278528, −507904, −147456, −180224), and shadow rows for weights 0–4 (88, 32, 18, 4, 2). These give exactly 32768(A6+A7) with right-hand side 65536 and zero residual in all 18 columns. The root session rechecked this identity directly against `split/system.py`.
- **Wording.** Infeasibility of the necessary relaxation excludes every code in that cell. Here, `cert=false` together with the message "infeasible" refers to the *dual* certificate search, so it does not show the primal is feasible. Sol solved the primal relaxations of the r=1 containment screen exactly and found them feasible.
- **local_exact.py.** The suffix-dimension row equals 4^9 × the local-rank row, so it adds nothing.
