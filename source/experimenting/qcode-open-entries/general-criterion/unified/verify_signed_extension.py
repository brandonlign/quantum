"""Compare the exact solver with independent enumeration of tiny Y matrices."""
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import random
from model import (ExtensionModel, coordinates, lift_coordinate, pullback_gram,
                   reduced_basis, recover_extension, span, embed, weight)
from signed_extension import solve_extension
from verify_unified import pauli

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]


def standard_columns(model, xs, ys):
    m, p = model.m, model.m - 1
    std = [1 << (m + i) for i in range(p)] + [(1 << i) | (1 << p) for i in range(p)]
    return [lift_coordinate(coordinates(e, model.R + std) >> 1, list(xs) + list(ys))
            for e in model.E]


def brute(model, xs):
    checked = 0
    for ys in itertools.product(range(1 << model.b), repeat=model.m - 1):
        checked += 1
        columns = standard_columns(model, xs, ys)
        if pullback_gram(columns, model.b) == model.Gamma and model.criterion(columns):
            return True, checked
    return False, checked


def positive_control():
    S = []
    for line in (ROOT / 'release16/search/q1645.txt').read_text().splitlines():
        x, z = line.split('|')
        S.append(sum(int(c) << i for i, c in enumerate(x.split())) |
                 (sum(int(c) << i for i, c in enumerate(z.split())) << 16))
    p = min((v for v in span(S) if v), key=lambda v: weight(v, 16))
    sites = [i for i in range(16) if ((p & 65535) | (p >> 16)) >> i & 1]
    rest = [i for i in range(16) if i not in sites]
    assert len(sites) == 6
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
        normalized.append(embed(px | (pz << 6), suffix, 6, 10))
    model, _ = recover_extension(normalized, 16, 4, 5, 6)
    leaf = json.loads((BASE / 'control_leaf.json').read_text())
    assert model.Gamma == leaf['Gamma'] and model.distances == leaf['dmin']
    result = solve_extension(model, leaf['xs'], node_limit=10000)
    assert result['status'] == 'FOUND' and model.criterion(result['model_columns'])
    return result


def verify():
    counts = Counter()
    matrices = 0
    cases = 0
    for generator in range(1, 16):
        for d in range(1, 6):
            model = ExtensionModel(5, 1, d, 3, [7 << 3], [generator])
            for xs in itertools.combinations(range(1, 8), 2):
                expected, checked = brute(model, xs)
                result = solve_extension(model, list(xs), node_limit=10000)
                assert result['status'] in ('FOUND', 'EXCLUDED')
                assert (result['status'] == 'FOUND') == expected
                if expected:
                    assert model.direct_property(result['model_columns'])
                counts[result['status']] += 1
                matrices += checked
                cases += 1
    rng = random.Random(1743)
    samples = []
    while len(samples) < 20:
        xs = rng.sample(range(1, 16), 3)
        if len(reduced_basis(xs)) == 3:
            samples.append(xs)
    for d in range(1, 6):
        model = ExtensionModel(6, 1, d, 4, [15 << 4], [])
        for xs in samples:
            expected, checked = brute(model, xs)
            result = solve_extension(model, xs, node_limit=10000)
            assert result['status'] in ('FOUND', 'EXCLUDED')
            assert (result['status'] == 'FOUND') == expected
            if expected:
                assert model.direct_property(result['model_columns'])
            counts[result['status']] += 1
            matrices += checked
            cases += 1
    small = ExtensionModel(5, 1, 2, 3, [7 << 3], [3 << 2])
    limited = solve_extension(small, [1, 2], node_limit=0)
    assert limited['status'] == 'LIMIT' and limited['nodes'] == 1
    control = positive_control()
    files = [Path(__file__), BASE / 'signed_extension.py', BASE / 'two_site.py',
             BASE / 'model.py', BASE / 'kernel_bound.py']
    report = dict(status='PASS', fixed_Z_cases=cases, brute_Y_matrices_checked=matrices,
                  solver_statuses=dict(counts), positive_control=control,
                  limit_control=limited, seed=1743,
                  sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                  scope='Brute Y enumeration uses no normal-form/gauge pruning; physical '
                        'checks share model.py. Small exact cases, known16,4,5 positive control '
                        'and LIMIT semantics; no complete new parameter exclusion.')
    (BASE / 'signed_extension_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('positive_control', 'sha256')}))


if __name__ == '__main__':
    verify()
