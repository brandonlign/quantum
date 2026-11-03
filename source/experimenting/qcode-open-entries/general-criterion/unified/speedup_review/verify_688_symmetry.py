"""Read-only check of the physical S6 action missing from the 688 Z search.

This is symmetry evidence, not an extension search or an exclusion certificate.
Run from unified with python3 -B speedup_review/verify_688_symmetry.py.
"""
import itertools
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model import coordinates, gram, lift_coordinate, reduced_basis, span, weight
from mixed_frames.frames import mixed_frame, xor_table


def permute(w, sites):
    return sum(((w >> j) & 1) << sites[j] for j in range(8)) | sum(
        ((w >> (8+j)) & 1) << (8+sites[j]) for j in range(8))


def normal(xs):
    # X0 is the tail singleton and stays fixed. The other five X images
    # label differences from the sixth equivalent support site.
    points = [0, *xs[1:]]
    return [xs[0], *min(sorted(v ^ a for v in points if v != a) for a in points)]


def main():
    f = mixed_frame('688')
    assert f.E[:6] == [0x4000, 0x2100, 0x2200, 0x2400, 0x2800, 0x3000]
    words = xor_table(f.E)
    physical_checks = 0
    # Adjacent swaps generate S6; the tail swap is an extra symmetry, but
    # it changes X0 and is intentionally omitted from the normal form.
    for j in range(5):
        sites = list(range(8))
        sites[j], sites[j+1] = sites[j+1], sites[j]
        assert set(permute(r, sites) for r in span(f.R)) == set(span(f.R))
        A = [coordinates(permute(w, sites), f.R + f.E) >> len(f.R) for w in f.E]
        assert len(reduced_basis(A)) == 12
        assert gram([permute(w, sites) for w in f.E], 8) == gram(f.E, 8)
        assert A[0] == 1
        assert all(a < 64 for a in A[:6])
        for e, w in enumerate(words):
            ee = lift_coordinate(e, A)
            assert f.costs[e] == f.costs[ee]
            assert sorted(weight(w ^ r, 8) for r in span(f.R)) == sorted(
                weight(words[ee] ^ r, 8) for r in span(f.R))
            physical_checks += 1
    actions = set()
    for order in itertools.permutations(range(6)):
        sites = [*order, 6, 7]
        A = tuple(coordinates(permute(w, sites), f.R + f.E) >> len(f.R) for w in f.E[:6])
        assert A[0] == 1 and all(a < 64 for a in A)
        actions.add(A)
    assert len(actions) == 720
    # Every independent map has exactly 720 distinct labeled images and
    # one anchor/sort normal form. Full Z searches enforce independence.
    fixture_checks = 0
    for xs in ([1, 2, 4, 8, 16, 32], [1370, 1, 2, 4, 8, 16]):
        assert len(reduced_basis(xs)) == 6
        orbit = {tuple(lift_coordinate(a, xs) for a in A) for A in actions}
        assert len(orbit) == 720
        assert {tuple(normal(x)) for x in orbit} == {tuple(normal(xs))}
        assert tuple(normal(xs)) in orbit
        fixture_checks += len(orbit)
    result = dict(status='PASS', scope='physical symmetry only; no search/exclusion',
                  quotient_checks=physical_checks, distinct_Z_actions=len(actions),
                  orbit_fixture_checks=fixture_checks,
                  tail_swap_used=False, runtime_speedup_measured=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
