"""Tiny fixture for NEXT.md: kernel-aware prefix capacity for local q=3 frames.

For each dimension-4 local H = S^perp cap V_T (containing Z^8) in
gauge_rankone_h5/all_critical_dim4.json, compute kappa(eta) = min weight on
coset eta of H in R^perp_T (2048 cosets, E/K of dimension 11) and the RHS
counts  #{eta != 0 : kappa(eta) >= 6-j}  of the THEOREM.md capacity
inequality.  Pauli encoding: integer g with x = g & 255, z = g >> 8.
"""
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
src = json.loads((ROOT / "gauge_rankone_h5/all_critical_dim4.json").read_text())
pc = lambda v: bin(v).count("1")
wt = lambda g: pc((g & 255) | (g >> 8))
sp = lambda a, b: (pc((a & 255) & (b >> 8)) + pc((a >> 8) & (b & 255))) & 1
Z8 = 255 << 8

def span(gens):
    out = [0]
    for g in gens:
        out += [w ^ g for w in out]
    return out

rows = []
for idx, cl in enumerate(src["classes"]):
    H = span(cl["generators"])
    assert len(set(H)) == 16
    if cl["gram_rank"] != 2:
        continue
    # census words are not normalized to Z^8; R=<p>, p the unique weight-8 word
    eights = [h for h in H if wt(h) == 8]
    assert len(eights) == 1
    p = eights[0]
    seen, kap = set(), Counter()
    for v in range(1 << 16):
        if sp(v, p) or v in seen:
            continue
        coset = [v ^ h for h in H]
        seen.update(coset)
        if v in H:
            continue
        kap[min(wt(c) for c in coset)] += 1
    rhs = {j: sum(n for k, n in kap.items() if k >= 6 - j) for j in range(1, 6)}
    rows.append({"class": idx, "gram_rank": cl["gram_rank"],
                 "spectrum": cl["spectrum"], "kappa_profile": dict(sorted(kap.items())),
                 "rhs_kappa_ge_6_minus_j": rhs})

for r in rows:
    print(r["class"], r["gram_rank"], r["kappa_profile"], r["rhs_kappa_ge_6_minus_j"])
(Path(__file__).with_suffix(".json")).write_text(json.dumps(rows, indent=1) + "\n")
