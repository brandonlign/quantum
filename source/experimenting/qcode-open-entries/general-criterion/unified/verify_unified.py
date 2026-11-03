"""Independent exact replay for the unified extension characterization.

Run with Python's standard library. No papers, release files or running
searches are modified. --output writes an optional JSON verification report.
"""
import argparse
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
from model import (ExtensionModel, coordinates, dual, embed, gram, lift_coordinate,
                   reduced_basis, recover_extension, span, split, symp, weight)

ROOT = Path(__file__).resolve().parents[2]


def pauli(text):
    n = len(text)
    return sum((c in 'XY') << i for i, c in enumerate(text)) | (
        sum((c in 'ZY') << i for i, c in enumerate(text)) << n)


def verify():
    report = {}
    # Enumerate ALL maps for ALL two-dimensional distance-two suffix codes.
    codes = {tuple(reduced_basis([u, v]))
             for u, v in combinations(range(1, 16), 2)}
    assert len(codes) == 35
    tested = accepted = weight_only_counterexamples = 0
    for C in sorted(codes):
        if min(weight(v, 2) for v in span(C) if v) < 2:
            continue
        model = ExtensionModel(4, 1, 2, 2, [pauli('ZZ')], list(C))
        for columns in product(range(1 << model.b), repeat=model.a):
            assert model.criterion(columns) == model.direct_property(columns)
            model.kernel_rank_identity(columns)
            tested += 1
            accepted += model.criterion(columns)
            weights = all(model.cost[e] + model.distances[model.map_value(e, columns)] >= 2
                          for e in range(1, 1 << model.a))
            if weights and not model.criterion(columns):
                weight_only_counterexamples += 1
    assert accepted and weight_only_counterexamples
    report['exhaustive_small_maps'] = dict(tested=tested, accepted=accepted,
        weight_admissible_incompatible=weight_only_counterexamples)

    # Rank-two local stabilizer: the known five-qubit [[5,1,3]] code.
    S5 = [pauli(s) for s in ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']]
    five, M5 = recover_extension(S5, 5, 1, 3, 4)
    assert five.r == 2 and five.criterion(M5) and five.direct_property(M5)
    q, t, g = five.kernel_rank_identity(M5)
    report['rank_two_pure_control'] = dict(parameters=[5, 1, 3], m=4, r=2,
        prefix_dimension=five.a, suffix_dimension=len(five.C),
        kernel_dimension=q, kernel_form_rank=t, gram_rank=g)

    # Localized impurity: append a fixed qubit to the same five-qubit code.
    S6 = [embed(pauli('Z'), 0, 1, 5)] + [embed(0, s, 1, 5) for s in S5]
    six, M6 = recover_extension(S6, 6, 1, 3, 1)
    assert six.r == 1 and six.a == 0 and six.criterion(M6)
    assert six.direct_property(M6) and min(weight(v, 6) for v in span(S6) if v) == 1
    report['localized_impure_control'] = dict(parameters=[6, 1, 3], m=1, r=1,
        minimum_stabilizer_weight=1, suffix_distance=six.C_distance)

    # Recover the complete 14-qubit obstruction from the general criterion.
    cells = json.loads((ROOT / 'release16/census/cells.json').read_text())
    reasons, detail = Counter(), []
    for i, C in enumerate(cells['10,10']):
        model = ExtensionModel(14, 3, 5, 4, [pauli('ZZZZ')], C)
        assert model.C_distance >= 5
        assert model.prefix_profile() == {0: 1, 1: 4, 2: 27, 3: 24, 4: 8}
        result = model.small_prefix_obstruction()
        assert result['excluded']
        reasons[result['reason']] += 1
        detail.append(dict(class_index=i, **result))
    assert dict(reasons) == {'hull': 25, 'capacity': 12}
    report['fourteen_qubit_replay'] = dict(classes=37, exclusions=dict(reasons),
        prefix_cost_profile={0: 1, 1: 4, 2: 27, 3: 24, 4: 8},
        capacities={1: 8, 2: 32, 3: 59}, per_class=detail)

    # Independently recover the existing [[16,4,5]] control without enumerating
    # its million-element centralizer. Normalize a minimum word to Z^6.
    S16 = []
    for line in (ROOT / 'release16/search/q1645.txt').read_text().splitlines():
        x, z = line.split('|')
        S16.append(sum(int(c) << i for i, c in enumerate(x.split())) |
                   (sum(int(c) << i for i, c in enumerate(z.split())) << 16))
    p = min((v for v in span(S16) if v), key=lambda v: weight(v, 16))
    T = [i for i in range(16) if ((p & 65535) | (p >> 16)) >> i & 1]
    assert len(T) == 6
    remaining = [i for i in range(16) if i not in T]
    normalized = []
    for v in S16:
        x, z = v & 65535, v >> 16
        px = pz = 0
        for j, i in enumerate(T):
            a, b = (p >> i) & 1, (p >> (16 + i)) & 1
            qx, qz = (1, 0) if b else (0, 1)
            u = (((x >> i) & 1) * b + ((z >> i) & 1) * a) & 1
            rem = (((x >> i) & 1) ^ (u * qx), ((z >> i) & 1) ^ (u * qz))
            assert rem in [(0, 0), (a, b)]
            px |= u << j
            pz |= int(rem != (0, 0)) << j
        w = sum(((x >> i) & 1) << j for j, i in enumerate(remaining)) | (
            sum(((z >> i) & 1) << j for j, i in enumerate(remaining)) << 10)
        normalized.append(embed(px | (pz << 6), w, 6, 10))
    control, Mc = recover_extension(normalized, 16, 4, 5, 6)
    assert control.r == 1 and control.criterion(Mc)
    control.kernel_rank_identity(Mc)
    # Serialize a standard Z_i/X_i X_last leaf for separate filter verification.
    Z = [1 << (6 + i) for i in range(5)]
    X = [(1 << i) | (1 << 5) for i in range(5)]
    xs = [control.map_value(coordinates(z, control.R + control.E) >> control.r, Mc)
          for z in Z]
    control_leaf = dict(xs=xs, dmin=control.distances, Gamma=control.Gamma,
                        b=control.b, m=6, Zbits=Z, Xbits=X)
    report['sixteen_qubit_control'] = dict(parameters=[16, 4, 5], m=6, r=1,
        syndrome_dimension=control.b, suffix_distance=control.C_distance,
        gram_rank=len(reduced_basis(control.Gamma)), exact_map_passes=True)

    report['m9_gram_ranks'] = dict(sorted(Counter(
        len(reduced_basis(gram(dual(C, 7), 7))) for C in cells['7,4']).items()))
    assert report['m9_gram_ranks'] == {6: 1, 8: 19, 10: 23}
    report['scope'] = ('Exact small-model converse checks, pure/rank-two/impure controls, '
                       'complete 14-qubit hull-capacity replay, and 16-qubit structural data. '
                       'This does not certify the pending m9 X-filter run or establish novelty.')
    return report, control_leaf


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--control-leaf', type=Path)
    args = parser.parse_args()
    report, leaf = verify()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    if args.control_leaf:
        args.control_leaf.write_text(json.dumps(leaf, indent=2) + '\n')
    for key, value in report.items():
        if key == 'fourteen_qubit_replay':
            value = {k: v for k, v in value.items() if k != 'per_class'}
        print(key + ':', value)
    print('UNIFIED VERIFICATION PASS')
