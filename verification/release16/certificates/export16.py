#!/usr/bin/env python3
"""Export [[16,5,5]] Farkas certificates (found by system.py) to the labelled JSON
format of the [[14,3,5]] paper's Appendix A, for checking by check_certificates.py."""
import json, sys
from fractions import Fraction
from math import lcm
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from system import build, farkas_exact, patterns

n, k, d = 16, 5, 5
CASES = [
    # name, blocks, weight equalities, weight >= (w,value) or None, meaning
    ('w1_stabilizer', [(1, ['Z'])], {}, None, 'S contains a weight-1 word (WLOG Z)'),
    ('w2_stabilizer', [(2, ['ZZ'])], {}, None, 'S contains a weight-2 word (WLOG ZZ)'),
    ('w3_stabilizer', [(3, ['ZZZ'])], {}, None, 'S contains a weight-3 word (WLOG ZZZ)'),
    ('w4_A4ge2', [(4, ['ZZZZ'])], {}, (4, 2), 'S contains ZZZZ and A4 >= 2'),
    ('w5_pure', [(5, ['ZZZZZ'])], {1: 0, 2: 0, 3: 0, 4: 0}, None,
     'A1=A2=A3=A4=0 and S contains a weight-5 word (WLOG ZZZZZ)'),
    ('minwt10', [], {w: 0 for w in range(1, 10)}, None, 'A1..A9 = 0 (stabilizer minimum weight >= 10)'),
]


def export(name, blocks, W, ge, meaning):
    weights = dict(W)
    if ge: weights[ge[0]] = ge[1]
    Meq, b, G, cols, le, ll = build(n, k, d, blocks, weights=weights or None)
    if ge:   # slack column turns the weight equality into  sum H - s = value, s >= 0
        Meq = np.hstack([Meq, np.zeros((Meq.shape[0], 1), dtype=object)]); Meq[-1, -1] = -1
        G = np.hstack([G, np.zeros((G.shape[0], 1), dtype=object)])
    cert, ok = farkas_exact(Meq, b, G)
    assert cert, name
    nu, lam = cert
    den = lcm(*[Fraction(v).denominator for v in list(nu) + list(lam) if v])
    rows = []
    def lab(L):
        if isinstance(L, str): return {'type': L}
        t, c = L
        if t == 'weight':
            w = c
            if ge and w == ge[0]: return None
            return {'type': 'weight', 'w': w, 'value': W[w]}
        return {'type': t, 'cell': list(c)}
    for i, v in enumerate(nu):
        if not v: continue
        m = int(v * den); r = lab(le[i])
        if r is None:   # the weight row carrying the >= assumption; nu must be <= 0
            assert m < 0
            r = {'type': 'weight_ge', 'w': ge[0], 'value': ge[1]}; m = -m
        r['multiplier'] = m; rows.append(r)
    for i, v in enumerate(lam):
        if not v: continue
        r = lab(ll[i]); r['multiplier'] = int(v * den); rows.append(r)
    # expected values (recomputed independently by the checker)
    res = [sum(int(nu[i] * den) * Meq[i, j] for i in range(len(nu)) if nu[i]) +
           sum(int(lam[i] * den) * G[i, j] for i in range(len(lam)) if lam[i])
           for j in range(len(cols))]
    nurhs = sum(int(nu[i] * den) * b[i] for i in range(len(nu)) if nu[i])
    blocks_json = []
    for m_, gens in blocks:
        pats = patterns(m_, gens)[1]
        blocks_json.append({'sites': m_, 'generators': gens, 'variables': 'coset_patterns',
                            'patterns': [[list(p), c] for p, c in sorted(pats.items())]})
    spec = {'name': name, 'code': '[[16,5,5]]', 'description': meaning,
            'blocks': blocks_json, 'suffix_sites': n - sum(m_ for m_, _ in blocks),
            'weight_assumptions': [[w, v] for w, v in W.items()],
            'weight_ge_assumptions': [list(ge)] if ge else [],
            'expected': {'columns': len(cols), 'nu_rhs': int(nurhs),
                         'positive_residuals': sum(1 for r in res if r > 0)},
            'certificate': rows}
    (HERE / f'{name}.json').write_text(json.dumps(spec, indent=1) + '\n')
    print(f'{name:16s} rows={len(rows):3d} columns={len(cols)} nu.rhs={nurhs}')


for c in CASES:
    export(*c)
