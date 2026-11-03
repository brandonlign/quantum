"""Distance-six positive control from Grassl's stored [[16,1,6]] matrix."""
import hashlib
import json
from pathlib import Path
from model import (coordinates, embed, reduced_basis, recover_extension, span, weight)
from signed_extension import solve_extension
from kernel_bound import kernel_screen

BASE = Path(__file__).resolve().parent
SOURCE = 'https://codetables.de/QECC.php?q=4&n=16&k=1'


def verify():
    data = BASE / 'positive_control_16_1_6.txt'
    S = []
    for line in data.read_text().splitlines():
        x, z = line.split('|')
        S.append(sum(int(c) << i for i, c in enumerate(x.split())) |
                 (sum(int(c) << i for i, c in enumerate(z.split())) << 16))
    assert len(reduced_basis(S)) == 15
    p = min((v for v in span(S) if v), key=lambda v: weight(v, 16))
    m = weight(p, 16)
    assert m == 6
    sites = [i for i in range(16) if ((p & 65535) | (p >> 16)) >> i & 1]
    rest = [i for i in range(16) if i not in sites]
    normalized = []
    for v in S:
        x, z = v & 65535, v >> 16
        px = pz = 0
        for j, i in enumerate(sites):
            a, b = (p >> i) & 1, (p >> (16 + i)) & 1
            qx, qz = (1, 0) if b else (0, 1)
            u = (((x >> i) & 1) * b + ((z >> i) & 1) * a) & 1
            rem = (((x >> i) & 1) ^ (u * qx), ((z >> i) & 1) ^ (u * qz))
            assert rem in [(0, 0), (a, b)]
            px |= u << j
            pz |= int(rem != (0, 0)) << j
        suffix = sum(((x >> i) & 1) << j for j, i in enumerate(rest)) | (
            sum(((z >> i) & 1) << j for j, i in enumerate(rest)) << 10)
        normalized.append(embed(px | (pz << m), suffix, m, 10))
    model, known_map = recover_extension(normalized, 16, 1, 6, m)
    assert model.r == 1 and model.criterion(known_map) and model.direct_property(known_map)
    xs = [model.map_value(coordinates(1 << (m + i), model.R + model.E) >> model.r,
                          known_map) for i in range(m - 1)]
    result = solve_extension(model, xs, node_limit=20000, minimum_stabilizer_weight=m)
    assert result['status'] == 'FOUND', result
    assert model.direct_property(result['model_columns'])
    found_stabilizer, _ = model.construct(result['model_columns'])
    actual_minimum = min(weight(v, 16) for v in span(found_stabilizer) if v)
    assert actual_minimum >= m
    bound = kernel_screen(16, 1, 6, m, 1, model.hull, minimum_stabilizer_weight=m)
    assert bound['status'] == 'UNRESOLVED'
    files = [data, Path(__file__), BASE / 'signed_extension.py', BASE / 'two_site.py',
             BASE / 'model.py', BASE / 'kernel_bound.py']
    report = dict(status='PASS', source=SOURCE, parameters=[16, 1, 6],
                  local_support=sites, m=m, r=model.r, suffix_dimension=len(model.C),
                  suffix_distance=model.C_distance, syndrome_dimension=model.b,
                  hull=model.hull, kernel_screen=bound, solver=result,
                  found_stabilizer_minimum=actual_minimum,
                  direct_centralizer_words_checked=1 << (model.n + model.k),
                  sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                  scope='Known positive control and reconstructed extension; not a new code '
                        'construction or a resolved open entry. Physical verifier shares model.py.')
    (BASE / 'distance6_control.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('solver','kernel_screen','sha256')}))


if __name__ == '__main__':
    verify()
