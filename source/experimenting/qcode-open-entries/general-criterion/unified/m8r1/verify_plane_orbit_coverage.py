"""Audit physical plane coverage and completed-record semantics, not search nodes."""
import hashlib
import itertools
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from suffix_symmetries import automorphisms

HERE = Path(__file__).resolve().parent


def physical_words(gens):
    words = {0}
    for g in gens:
        words |= {w ^ g for w in tuple(words)}
    return words


def weight(w):
    return ((w & 511) | (w >> 9)).bit_count()


def main():
    checkpoint = HERE / 'multi_plane_pilot.json'
    data = json.loads(checkpoint.read_text())
    for path, expected in data['sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, path
    classes = json.loads((HERE.parent / 'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
    plan = json.loads((HERE / 'heavy_plane_orbits.json').read_text())['records']
    closed = {int(i) for i, status in data['statuses'].items() if status == 'EXCLUDED'}
    checked_planes = checked_orbits = leaves = 0
    for item in plan:
        index = item['cls']
        if index not in closed:
            continue
        C = classes[index]['representative']['generators']
        heavy = {w for w in physical_words(C) if weight(w) >= 8}
        # Independently enumerate every two-dimensional physical plane.
        planes = {tuple(sorted((a, b, a ^ b))) for a, b in itertools.combinations(heavy, 2)
                  if a ^ b in heavy}
        assert len(planes) == item['plane_count']
        maps, _ = automorphisms(C, 9)
        remaining = set(planes)
        for orbit_index, orbit in enumerate(item['orbits']):
            a, b = orbit['representative']
            seed = tuple(sorted((a, b, a ^ b)))
            actual, frontier = {seed}, [seed]
            while frontier:
                plane = frontier.pop()
                for transform in maps:
                    image = tuple(sorted(transform(w) for w in plane))
                    assert image in planes
                    if image not in actual:
                        actual.add(image)
                        frontier.append(image)
            assert len(actual) == orbit['size']
            assert actual <= remaining, 'Overlapping/missing physical plane orbit'
            remaining -= actual
            record = data['records'][f'{index}:{orbit_index}']
            assert tuple(sorted((record['N'][0], record['N'][1], record['N'][0] ^ record['N'][1]))) == seed
            status = record['status']
            if status == 'EXCLUDED':
                assert record['z']['status'] == 'COMPLETE'
                assert not any(count for verdict, count in record['tally'].items()
                               if verdict.split(':', 1)[0] != 'EXCLUDED')
                assert record['found'] is None
                assert sum(record['tally'].values()) == record['z']['leaves']
                leaves += record['z']['leaves']
            elif status == 'EXCLUDED_RANK_TWO_PROPAGATION':
                supports = [(w & 511) | (w >> 9) for w in seed]
                assert supports[0] == supports[1] == supports[2]
                assert supports[0].bit_count() == 8
                r2 = json.loads((HERE.parent / 'rank2_branch_result.json').read_text())
                assert r2['status'] == 'COMPLETE_BRANCH_EXCLUSION' and r2['local_branch'] == {'m': 8, 'r': 2}
            elif status == 'EXCLUDED_FULL_SUPPORT_AUDIT':
                audited = json.loads((HERE / 'full_support_geometry_audit.json').read_text())['records']
                assert any(r['cls'] == index and r['status'] == 'EXCLUDED_ALTERNATE' for r in audited)
            else:
                raise AssertionError(f'Incomplete record counted as exclusion: {status}')
            checked_orbits += 1
        assert not remaining, 'Uncovered physical planes'
        checked_planes += len(planes)
    result = dict(status='PASS', closed_classes=len(closed), physical_planes=checked_planes,
                  completed_orbits=checked_orbits, checked_leaf_verdicts=leaves,
                  checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                  scope='Independent physical-plane enumeration and orbit BFS; shared verified automorphism generator; metadata audit, not exhaustive node replay')
    (HERE / 'plane_orbit_coverage_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(result)


if __name__ == '__main__':
    main()
