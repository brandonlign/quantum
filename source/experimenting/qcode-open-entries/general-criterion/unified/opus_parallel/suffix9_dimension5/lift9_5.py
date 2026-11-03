"""Complete one-coordinate lift (8,2^5,>=5) -> (9,2^5,>=6), then target17 m8,r1 filters.

Parents: cell_8_5.json (26299 classes, X|Z layout x|z<<8). For each parent, all
187 canonical tail-form pairs (1 rank0, 31 rank1, 155 rank2; PARENT_LIFT.md);
valid iff every weight-5 parent coefficient word gets a nonzero new symbol.
Then: hull, radical weight profile, maximal dimension of a radical subspace with
all nonzero weights >=8 (needed for N=ker tau <= rad C, min weight 8), and
exact nauty equivalence classes (site perms + local GL(2,2)).
Checkpoints every 2000 parents to lift9_5_checkpoint.json.
"""
import json, time, hashlib
from collections import Counter
from itertools import combinations
from pathlib import Path
import numpy as np
from pynauty import Graph, certificate

OUT = Path(__file__).resolve().parent
L8 = (1 << 8) - 1


def wt(w, n): return ((w & ((1 << n) - 1)) | (w >> n)).bit_count()
def sp(u, v, n):
    m = (1 << n) - 1
    return (((u & m) & (v >> n)).bit_count() + ((v & m) & (u >> n)).bit_count()) & 1
def rank(rows):
    piv = {}
    for r in rows:
        while r:
            b = r.bit_length() - 1
            if b in piv: r ^= piv[b]
            else: piv[b] = r; break
    return len(piv)
def words(gens):
    out = [0]
    for g in gens: out += [w ^ g for w in out]
    return out


def maps(c):
    out = [(0, 0)] + [(a, 0) for a in range(1, 1 << c)]
    out += [(a, b) for a in range(1, 1 << c) for b in range(a + 1, 1 << c) if b < a ^ b]
    return out


def lift(gens, a, b, f=8):
    m = (1 << f) - 1
    return [(g & m) | (((a >> i) & 1) << f) | (((g >> f) | (((b >> i) & 1) << f)) << (f + 1)) for i, g in enumerate(gens)]


def cert(gens, n):
    ws = words(gens); cnt = len(ws); adj = {i: set() for i in range(cnt + 3 * n)}
    for row, w in enumerate(ws):
        for j in range(n):
            s = ((w >> j) & 1) | (((w >> (n + j)) & 1) << 1)
            if s: adj[row].add(cnt + 3 * j + s - 1)
    for j in range(n):
        for s in range(3): adj[cnt + 3 * j + s].update(cnt + 3 * j + t for t in range(3) if t != s)
    return certificate(Graph(cnt + 3 * n, adjacency_dict=adj,
                             vertex_coloring=[set(range(cnt)), set(range(cnt, cnt + 3 * n))])).hex()


def heavy_dim(rad, n, thr=8):
    """Largest dim of a subspace of span(rad) whose nonzero words all have weight>=thr."""
    heavy = [w for w in rad if w and wt(w, n) >= thr]; best = 0
    for d in range(1, rank(rad) + 1):
        ok = False
        for basis in combinations(heavy, d):
            if rank(list(basis)) == d and all(wt(w, n) >= thr for w in words(list(basis)) if w):
                ok = True; break
        if ok: best = d
        else: break
    return best


def main():
    t0 = time.monotonic(); cell = json.loads((OUT / 'cell_8_5.json').read_text())
    assert cell['status'] == 'COMPLETE' and cell['mass_equal']
    parents = cell['xz']; forms = maps(5)
    assert len(forms) == 187 and sum(1 for a, b in forms if b) == 155
    par = np.array([[(a & w).bit_count() & 1 for w in range(32)] for a in range(32)], dtype=bool)
    cov = np.asarray([par[a] | par[b] for a, b in forms])
    lifts = []; per = Counter()
    for i, gens in enumerate(parents):
        ws = words(gens); assert len(set(ws)) == 32 and min(wt(w, 8) for w in ws if w) >= 5
        bad = [j for j, w in enumerate(ws) if wt(w, 8) == 5]
        for idx in np.flatnonzero(cov[:, bad].all(axis=1)):
            a, b = forms[int(idx)]; new = lift(gens, a, b)
            nw = words(new); assert len(set(nw)) == 32
            dist = min(wt(w, 9) for w in nw if w); assert dist >= 6
            G = [sum(sp(u, v, 9) << j for j, v in enumerate(new)) for u in new]
            h = 5 - rank(G)
            rad = [w for w in nw if all(sp(w, g, 9) == 0 for g in new)]; assert len(rad) == 1 << h
            lifts.append(dict(parent=i, tail_forms=[a, b], generators=new, distance=dist, hull=h,
                              radical_weights=sorted(wt(w, 9) for w in rad if w),
                              heavy8_radical_dim=heavy_dim(rad, 9),
                              weight_distribution=dict(sorted(Counter(wt(w, 9) for w in nw).items()))))
        if i % 2000 == 0:
            (OUT / 'lift9_5_checkpoint.json').write_text(json.dumps(dict(status='RUNNING', parents_done=i + 1, lifts=len(lifts))) + '\n')
    classes = {}
    for j, it in enumerate(lifts):
        classes.setdefault(cert(it['generators'], 9), []).append(j)
    cls = [dict(class_index=c, certificate_hex=k, representative=lifts[v[0]], extension_indices=v)
           for c, (k, v) in enumerate(sorted(classes.items()))]
    for c in cls:  # invariants must agree inside a class
        for j in c['extension_indices']:
            assert lifts[j]['hull'] == c['representative']['hull']
            assert lifts[j]['radical_weights'] == c['representative']['radical_weights']
    # m8,r1 (KERNEL_BOUND.md / kernel_bound_verification.json): (q,t,hull,dimN)
    poss = [(1, 0, 1, 0), (2, 0, 3, 1), (2, 2, 1, 1), (3, 0, 5, 2), (3, 2, 3, 2)]
    target = []
    for c in cls:
        r = c['representative']
        ok = [dict(q=q, t=t, hull=h, dimN=n) for q, t, h, n in poss if r['hull'] == h and r['heavy8_radical_dim'] >= n]
        if ok: target.append(dict(class_index=c['class_index'], hull=r['hull'], radical_weights=r['radical_weights'],
                                  heavy8_radical_dim=r['heavy8_radical_dim'], allowed_qtN=ok))
    res = dict(status='COMPLETE', parent_cell='8,5', parent_count=len(parents), canonical_forms_per_parent=187,
               extension_representatives=len(lifts), hull_counts_reps=dict(Counter(l['hull'] for l in lifts)),
               inequivalent_classes=len(cls), hull_counts_classes=dict(Counter(c['representative']['hull'] for c in cls)),
               target_m8r1_possibilities=[dict(q=q, t=t, hull=h, dimN=n) for q, t, h, n in poss],
               target_surviving_classes=len(target),
               target_surviving_by_qtN=dict(Counter(f"q{o['q']}t{o['t']}N{o['dimN']}" for x in target for o in x['allowed_qtN'])),
               target=target, classes=cls, seconds=round(time.monotonic() - t0, 1),
               sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT / 'cell_8_5.json', Path(__file__)]},
               scope='Complete relative to cell_8_5.json. Necessary filters only; no quantum extension tested.')
    (OUT / 'suffix9_dim5.json').write_text(json.dumps(res) + '\n')
    (OUT / 'lift9_5_checkpoint.json').write_text(json.dumps(dict(status='COMPLETE', lifts=len(lifts))) + '\n')
    print(json.dumps({k: v for k, v in res.items() if k not in ['classes', 'target', 'sha256']}))


if __name__ == '__main__':
    main()
