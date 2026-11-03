"""Unrestricted replay of the lift: all 1024 ordered tail-form pairs per (8,2^5) parent,
orbit multiplicities 1/3/6, regenerated generators, distance, hull, radical weights."""
import json
from collections import Counter
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parent


def pw(w, n): return ((w & ((1 << n) - 1)) | (w >> n)).bit_count()
def sp(a, b, n):
    m = (1 << n) - 1
    return (((a & m) & (b >> n)).bit_count() + ((b & m) & (a >> n)).bit_count()) & 1
def rank(rows):
    p = {}
    for r in rows:
        while r:
            b = r.bit_length() - 1
            if b in p: r ^= p[b]
            else: p[b] = r; break
    return len(p)
def canon(a, b):
    v = sorted({a, b, a ^ b} - {0})
    return (0, 0) if not v else (v[0], 0) if len(v) == 1 else tuple(v[:2])


cell = json.loads((OUT / 'cell_8_5.json').read_text())['xz']
res = json.loads((OUT / 'suffix9_dim5.json').read_text())
saved = {}
for c in res['classes']:
    pass
lifts = {}
for c in res['classes']:
    for j in c['extension_indices']: lifts[j] = None
# rebuild full lift list order from classes is not stored; regenerate canonical set per parent instead
forms = [(a, b) for a in range(32) for b in range(32)]
par = np.array([[(a & w).bit_count() & 1 for w in range(32)] for a in range(32)], dtype=bool)
cov = np.asarray([par[a] | par[b] for a, b in forms])
total = 0; unrestricted = 0; hulls = Counter()
for i, gens in enumerate(cell):
    ws = [0]
    for g in gens: ws += [w ^ g for w in ws]
    bad = [j for j, w in enumerate(ws) if pw(w, 8) == 5]
    good = np.flatnonzero(cov[:, bad].all(axis=1)); unrestricted += len(good)
    orb = Counter(canon(*forms[int(x)]) for x in good)
    for (a, b), mult in orb.items():
        assert mult == (1 if a == b == 0 else 3 if b == 0 else 6)
        new = [(g & 255) | (((a >> k) & 1) << 8) | ((g >> 8) << 9) | (((b >> k) & 1) << 17) for k, g in enumerate(gens)]
        e = [0]
        for g in new: e += [w ^ g for w in e]
        assert len(set(e)) == 32 and min(pw(w, 9) for w in e if w) >= 6
        G = [sum(sp(u, v, 9) << j for j, v in enumerate(new)) for u in new]
        hulls[5 - rank(G)] += 1; total += 1
out = dict(canonical_extensions=total, unrestricted_ordered_pairs=unrestricted, hull_counts=dict(hulls),
           matches_lift=total == res['extension_representatives'] and dict(hulls) == {int(k): v for k, v in res['hull_counts_reps'].items()})
(OUT / 'verify_lift9_5.json').write_text(json.dumps(out, indent=2) + '\n'); print(out)
assert out['matches_lift']
