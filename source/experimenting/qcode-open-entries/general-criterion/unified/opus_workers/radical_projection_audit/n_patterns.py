"""Tiny fixture for NEXT.md.

(1) Which census rank-2 dim-4 local planes have their weight-8 word p in rad(H)?
    (Symplectic product and weight are invariant under swapping x/z halves,
    so the test does not depend on the bit convention.)
(2) Support patterns of a 2-dim isotropic additive code N=<a,b> on the 9
    complement sites (q=3 case: dim S cap V_W = q-1 = 2), min weight >= 6.
    Per site: n0 zero, n2 sites of local dim 2 (all three words nonzero),
    and dim-1 sites where exactly one of a, b, a+b vanishes (na, nb, nc).
    Isotropy forces n2 even. Rules applied:
      R1  no word has weight < 6;
      R2  at most two words of weight 6 or 7 (A6+A7 <= 2);
      R3  a weight-8 word's support contains no other nonzero word's support
          (exclusive rank-one column: local rank on that support is 1).
"""
import itertools, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
pc = lambda v: bin(v).count("1")
wt = lambda g: pc((g & 255) | (g >> 8))
sp = lambda a, b: (pc(a & 255 & (b >> 8)) + pc((a >> 8) & b & 255)) & 1

out = {"census_rank2": [], "n_patterns": []}
src = json.loads((ROOT / "gauge_rankone_h5/all_critical_dim4.json").read_text())
for i, c in enumerate(src["classes"]):
    if c["gram_rank"] != 2:
        continue
    H = [0]
    for g in c["generators"]:
        H += [w ^ g for w in H]
    p = [h for h in H if wt(h) == 8]
    radical = all(sp(p[0], h) == 0 for h in H)
    xw = None
    if radical:
        # X-part relative to p: site letters equal to p's letter or I count as non-X
        def xbits(h):
            v = 0
            for s in range(8):
                hs = ((h >> s) & 1, (h >> (8 + s)) & 1)
                ps = ((p[0] >> s) & 1, (p[0] >> (8 + s)) & 1)
                if hs != (0, 0) and hs != ps:
                    v |= 1 << s
            return v
        xw = sorted({(xbits(h), wt(h)) for h in H if h not in (0, p[0])})
        xw = sorted(pc(x) for x, _ in {(xbits(h) , 0) for h in H if h not in (0, p[0])})
    out["census_rank2"].append({"class": i, "p_in_radical": radical, "x_weights_of_K": xw})

for n2 in range(0, 10, 2):
    for na, nb, nc in itertools.product(range(10), repeat=3):
        n0 = 9 - n2 - na - nb - nc
        if n0 < 0:
            continue
        # supports as site-type sets; a vanishes on na-sites, b on nb, a+b on nc
        S = {"a": {"2", "b", "c"}, "b": {"2", "a", "c"}, "c": {"2", "a", "b"}}
        cnt = {"2": n2, "a": na, "b": nb, "c": nc}
        w = {k: sum(cnt[t] for t in v) for k, v in S.items()}
        if min(w.values()) < 6:
            continue
        light = sum(1 for v in w.values() if v <= 7)
        if light > 2:
            continue
        def contained(x, y):  # supp x subset of supp y
            return all(cnt[t] == 0 for t in S[x] - S[y])
        bad = any(w[y] == 8 and contained(x, y) for x in S for y in S if x != y)
        if bad:
            continue
        key = tuple(sorted(w.values()))
        out["n_patterns"].append({"n0": n0, "n2": n2, "na": na, "nb": nb, "nc": nc, "weights": key})

pats = sorted({tuple(p["weights"]) for p in out["n_patterns"]})
out["weight_triples"] = pats
print(json.dumps(out["census_rank2"]))
print(pats, len(out["n_patterns"]))
Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=1) + "\n")
