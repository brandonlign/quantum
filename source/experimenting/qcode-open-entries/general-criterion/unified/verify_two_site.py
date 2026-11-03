"""Positive control, physical prefix costs and sampled congruence checks."""
import json
from pathlib import Path
import random
from two_site import screen, requirements, normal_form
from model import coordinates, reduced_basis

BASE = Path(__file__).resolve().parent


def verify():
    ctrl = json.loads((BASE / 'control_leaf.json').read_text())
    result = screen(ctrl['xs'], ctrl['dmin'], ctrl['Gamma'], ctrl['m'])
    assert result['verdict'] == 'PASS', result
    count = 0
    for m in range(3, 11):
        for support in range(1 << m):
            if support.bit_count() % 2:
                continue
            for d in [3, 4, 5, 6]:
                assert requirements(m, support, d) == requirements(m, support, d, True)
                count += 1
    rng = random.Random(4132)
    forms = 0
    for b in range(2, 8):
        for a in range(1, b + 1):
            for trial in range(10):
                X = []
                while len(X) < a:
                    value = rng.randrange(1, 1 << b)
                    if len(reduced_basis(X + [value])) > len(X):
                        X.append(value)
                Y = [rng.randrange(1 << b) for _ in range(a)]
                G = [0] * b
                for x, y in zip(X, Y):
                    for i in range(b):
                        G[i] ^= (y if x >> i & 1 else 0) ^ (x if y >> i & 1 else 0)
                nf = normal_form(X, G)
                yp = [coordinates(y, nf['T']) for y in Y]
                mask = (1 << a) - 1
                assert all(y & ~mask == nf['offsets'][i] for i, y in enumerate(yp))
                assert all(((((yp[i] >> j) ^ (yp[j] >> i)) & 1)
                            == ((nf['Gamma'][i] >> j) & 1))
                           for i in range(a) for j in range(a))
                assert len(reduced_basis(X + Y)) == a + len(reduced_basis(
                    [offset >> a for offset in nf['offsets']]))
                forms += 1
    return dict(positive_control=result['verdict'], physical_cost_checks=count,
                normal_form_checks=forms, status='PASS', seed=4132)


if __name__ == '__main__':
    report = verify()
    (BASE / 'two_site_controls.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))
