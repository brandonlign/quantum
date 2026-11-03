"""Compare relative profiles with direct physical logical distance."""
from collections import Counter
import hashlib
import itertools
import json
from math import inf
from pathlib import Path
from model import ExtensionModel, recover_extension, span, weight, symp
from relative_extension import logical_distance
from verify_unified import pauli

BASE = Path(__file__).resolve().parent


def direct(model, columns):
    S, D = model.construct(columns)
    assert all(symp(u, v, model.n) == 0 for u in S for v in S)
    stabilizer = set(span(S))
    return min((weight(v, model.n) for v in span(D) if v not in stabilizer), default=inf)


def verify():
    total = compatible = 0
    distances = Counter()
    for generator in range(1, 16):
        model = ExtensionModel(5, 1, 1, 3, [7 << 3], [generator])
        for columns in itertools.product(range(8), repeat=4):
            total += 1
            result = logical_distance(model, columns)
            if not result['compatible']:
                continue
            compatible += 1
            actual = direct(model, columns)
            assert result['distance'] == actual
            distances[actual] += 1
    S5 = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
    five, five_map = recover_extension(list(map(pauli,S5)),5,1,3,4)
    five_profile = logical_distance(five,five_map)
    assert five.r==2 and five_profile['distance']==direct(five,five_map)==3
    S7 = [pauli('ZIIIIII'), pauli('IZIIIII')] + [pauli('II'+s) for s in S5]
    impure, columns = recover_extension(S7, 7, 1, 3, 1)
    result = logical_distance(impure, columns)
    assert impure.C_distance == 1
    assert not impure.criterion(columns)
    assert result['compatible'] and result['distance'] == direct(impure, columns) == 3
    files = [Path(__file__), BASE / 'relative_extension.py', BASE / 'model.py']
    report = dict(status='PASS', maps_examined=total, compatible_maps_checked=compatible,
                  logical_distances=dict(distances),
                  rank_two_positive_control=dict(parameters=[5,1,3],m=4,r=2,distance=3),
                  nonlocalized_impure_control=dict(parameters=[7,1,3], m=1,r=1,
                        suffix_distance=impure.C_distance, logical_distance=result['distance'],
                        localized_criterion_passes=False, relative_criterion_passes=True),
                  sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                  scope='All maps in the tiny models, with suffix profile vs direct physical '
                        'logical-distance comparison. Construction and algebra share model.py; '
                        'no novelty or large impure-case exclusion claim.')
    (BASE / 'relative_extension_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))


if __name__ == '__main__':
    verify()
