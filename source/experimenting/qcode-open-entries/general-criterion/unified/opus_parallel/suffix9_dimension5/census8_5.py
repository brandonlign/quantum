"""Complete (8,2^5,>=5)_4 cell from the complete, mass-checked parent cells
(7,2^5) [r=0], (7,2^4) [r=1], (7,2^3) [r=2] produced by census85.py.

Extension rule identical to release16/census/mass_census.py (exec'd read-only):
new column rank r, r-dim subspace U of F2^14/C with every nonzero element of
coset-minimum weight >=4. Labelled extension counts ext_r are computed from the
full deep-coset list (no reduction). To limit nauty calls, U is reduced to one
representative per Aut(C)-orbit (Aut(C) generators from nauty, acting on the
quotient); this reduction is validated a posteriori by the labelled mass
identity  sum_P (7!6^7/|Aut P|) ext_r(P) = sum_D 8!6^8/|Aut D|.

Interleaved layout internally; output also in release16 X|Z layout.
Usage: python3.11 census8_5.py count | run
"""
import json, sys, time
from fractions import Fraction
from math import factorial
from pathlib import Path
from pynauty import autgrp

OUT = Path(__file__).resolve().parent
SRC = Path('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/mass_census.py')
ns = {'__file__': str(SRC)}
exec(SRC.read_text().split('if __name__ == "__main__":')[0], ns)
span, wt, graph, aut_order, certificate = (ns[k] for k in ['span', 'wt', 'graph', 'aut_order', 'certificate'])
quotient_labels, light_errors = ns['quotient_labels'], ns['light_errors']
from census85 import to_xz


def parents():
    out = {}
    for f, key in [('cells_d5_7_5.json', '7,5'), ('cells_d5_7_4.json', '7,4'), ('cells_d5_7_3.json', '7,3')]:
        d = json.loads((OUT / f).read_text())[key]
        assert d['mass_equal'] and len(d['interleaved_with_aut']) == d['classes']
        out[5 - int(key[2])] = d['interleaved_with_aut']
    return out


def deep_cosets(gens, n):
    label, lift, q = quotient_labels(gens, n)
    covered = {label(e) for e in light_errors(n, 3)}
    return label, lift, [a for a in range(1 << q) if a not in covered]


def monomials(gens, n):
    W = len(span(gens)); out = []
    for p in autgrp(graph(gens, n))[0]:
        pi = [0] * n; sig = [[0] * 4 for _ in range(n)]
        for i in range(n):
            for s in (1, 2, 3):
                t = p[W + 3 * i + s - 1] - W
                pi[i] = t // 3; sig[i][s] = t % 3 + 1
        out.append((pi, sig))
    return out


def apply(m, v, n):
    pi, sig = m; out = 0
    for i in range(n):
        s = (v >> (2 * i)) & 3
        if s: out |= sig[i][s] << (2 * pi[i])
    return out


def find(par, x):
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x


def orbit_reps(items, images):
    par = list(range(len(items)))
    for img in images:
        for j, y in enumerate(img):
            a, b = find(par, j), find(par, y)
            if a != b: par[a] = b
    return sorted({find(par, j) for j in range(len(items))})


def extend(gens, r, n=7):
    label, lift, deep = deep_cosets(gens, n)
    ds = set(deep); top = 2 * n
    if r == 0:
        return 1, [list(gens)]
    monos = monomials(gens, n)
    act = [{s: label(apply(m, lift(s), n)) for s in deep} for m in monos]
    if r == 1:
        idx = {s: j for j, s in enumerate(deep)}
        reps = orbit_reps(deep, [[idx[a[s]] for s in deep] for a in act])
        return 3 * len(deep), [list(gens) + [lift(deep[j]) | (1 << top)] for j in reps]
    subs = []
    for i, a in enumerate(deep):
        for b in deep[i + 1:]:
            x = a ^ b
            if x > b and x in ds: subs.append((a, b, x))
    idx = {s: j for j, s in enumerate(subs)}
    imgs = [[idx[tuple(sorted((a[u], a[v], a[w])))] for (u, v, w) in subs] for a in act]
    reps = orbit_reps(subs, imgs)
    return 6 * len(subs), [list(gens) + [lift(subs[j][0]) | (1 << top), lift(subs[j][1]) | (2 << top)] for j in reps]


def count():
    for r, ps in sorted(parents().items()):
        tot = 0
        for gens, aut in ps:
            _, _, deep = deep_cosets(gens, 7); ds = set(deep)
            sub = sum(1 for i, a in enumerate(deep) for b in deep[i + 1:] if (a ^ b) > b and (a ^ b) in ds)
            tot += [1, len(deep), sub][r]
        print(json.dumps(dict(r=r, parents=len(ps), unreduced_subspaces=tot)), flush=True)


def run():
    t0 = time.time(); ps = parents(); classes = {}; lhs = Fraction(0); reps_total = 0; per = []
    for r in (0, 1, 2):
        for i, (gens, aut) in enumerate(ps[r]):
            ext, gs = extend(gens, r)
            lhs += Fraction(factorial(7) * 6 ** 7, aut) * ext
            reps_total += len(gs)
            per.append(dict(r=r, parent=i, labelled_ext=ext, orbit_reps=len(gs)))
            for g in gs:
                G = graph(g, 8); c = certificate(G)
                if c not in classes: classes[c] = (g, aut_order(G), r, i)
    found = list(classes.values())
    for g, a, _, _ in found:
        sp = span(g); assert len(set(sp)) == 32 and min(wt(v, 8) for v in sp if v) >= 5
        assert (factorial(8) * 6 ** 8) % a == 0
    rhs = sum(Fraction(factorial(8) * 6 ** 8, a) for g, a, _, _ in found)
    res = dict(status='COMPLETE' if lhs == rhs else 'MASS_MISMATCH', cell='8,5', classes=len(found),
               orbit_representatives_tested=reps_total, labelled_mass_lhs=str(lhs), labelled_mass_rhs=str(rhs),
               mass_equal=lhs == rhs, seconds=round(time.time() - t0, 1), per_parent=per,
               interleaved=[g for g, *_ in found], aut_orders=[a for _, a, _, _ in found],
               origin=[[r, i] for *_, r, i in found], xz=[[to_xz(v, 8) for v in g] for g, *_ in found],
               layout_xz='word = x | z<<8 (release16 cells.json convention)')
    (OUT / 'cell_8_5.json').write_text(json.dumps(res) + '\n')
    print(json.dumps({k: v for k, v in res.items() if k not in ['per_parent', 'interleaved', 'xz', 'aut_orders', 'origin']}))
    assert lhs == rhs


if __name__ == '__main__':
    {'count': count, 'run': run}[sys.argv[1]]()
