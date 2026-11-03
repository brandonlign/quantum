"""Remove release-census assumptions using exact parent/child incidences."""
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from suffix_mass_check.verify8 import type_check
from suffix_mass_check.common import canonical, coset_minima
from model import reduced_basis, span, weight
from mixed_frames.stage3.common import atomic_json, sha


def main():
    start = time.monotonic()
    data = ROOT / 'opus_parallel/suffix9_dimension5'
    paths = [data / f'cells_d5_7_{k}.json' for k in (3, 4, 5)]
    frozen = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    cells = {k: json.loads(p.read_text())[f'7,{k}']['xz'] for k, p in zip((3, 4, 5), paths)}
    assert [len(cells[k]) for k in (3, 4, 5)] == [32, 43, 1]
    types = type_check(cells[3])
    assert types == dict(multisets=116280, feasible=2171, classes=32, star_certificate_sets_equal=True)
    print('Independent dimension-three coordinate-type enumeration PASS', flush=True)
    group = math.factorial(7) * 6 ** 7
    info = {}
    for k in (3, 4, 5):
        seen = set()
        info[k] = []
        for index, C in enumerate(cells[k]):
            assert len(reduced_basis(C)) == k and min(weight(w, 7) for w in span(C) if w) >= 5
            key, order = canonical(C, 7, need_order=True)
            assert key not in seen and group % order == 0
            seen.add(key)
            distances = coset_minima(C, 7)
            info[k].append(dict(cls=index, certificate=key, orbit_size=group // order,
                                automorphisms=order, extensions=int((distances >= 5).sum())))
    transitions = []
    for k in (3, 4):
        lookup = {r['certificate']: r['cls'] for r in info[k]}
        weighted = Counter()
        hyperplanes = 0
        for D, record in zip(cells[k + 1], info[k + 1]):
            for functional in range(1, 1 << (k + 1)):
                pivot = (functional & -functional).bit_length() - 1
                C = [D[j] ^ (D[pivot] if functional & (1 << j) else 0)
                     for j in range(k + 1) if j != pivot]
                assert len(reduced_basis(C)) == k
                key = canonical(C, 7)
                assert key in lookup
                weighted[lookup[key]] += record['orbit_size']
                hyperplanes += 1
        for record in info[k]:
            assert weighted[record['cls']] == record['orbit_size'] * record['extensions'], (k, record['cls'])
        lhs = sum(r['orbit_size'] * r['extensions'] for r in info[k])
        mass = sum(r['orbit_size'] for r in info[k + 1])
        assert lhs == mass * ((1 << (k + 1)) - 1)
        transitions.append(dict(parent_dimension=k, child_dimension=k + 1,
                                parent_classes=len(cells[k]), child_classes=len(cells[k + 1]),
                                incidence=lhs, child_physical_mass=mass,
                                every_parent_pointwise_incidence_matches=True,
                                hyperplanes_enumerated=hyperplanes))
    assert frozen == {str(p.relative_to(ROOT)): sha(p) for p in paths}
    records = {str(k): [{name: value for name, value in r.items() if name != 'certificate'}
                        for r in info[k]] for k in (3, 4, 5)}
    report = dict(status='PASS_COMPLETE_7_SITE_D5_PARENTS_FROM_TYPES',
                  dimension3_type_enumeration=types, transitions=transitions, records=records,
                  seconds=time.monotonic() - start,
                  hashes=dict(frozen, **{str(p.relative_to(ROOT)): sha(p) for p in
                              [Path(__file__), HERE / 'PROOF.md', ROOT / 'suffix_mass_check/verify8.py',
                               ROOT / 'suffix_mass_check/common.py', ROOT / 'model.py', ROOT / 'mixed_frames/frames.py']}),
                  scope='Complete 7-site dimensions 3,4,5 distance-at-least-five additive-code cells from exhaustive coordinate types and exact pointwise incidence; no release-census completeness assumption, no quantum extension exclusion')
    atomic_json(HERE / 'verification.json', report)
    print(json.dumps({k: v for k, v in report.items() if k not in ('records', 'hashes')}), flush=True)


if __name__ == '__main__':
    main()
