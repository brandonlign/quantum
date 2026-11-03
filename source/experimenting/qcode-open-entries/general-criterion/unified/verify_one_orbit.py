"""Nonvacuous small extension control for the one-orbit corollary."""
import itertools
import json
from pathlib import Path
from model import ExtensionModel, coordinates, lift_coordinate
from two_site import normal_form, screen


def verify():
    m = 3
    model = ExtensionModel(5, 1, 2, m, [7 << m], [3 << 2])
    xs = [1, 2]
    assert screen(xs, model.distances, model.Gamma, m, d=2)['verdict'] == 'PASS'
    nf = normal_form(xs, model.Gamma)
    std = [1 << m, 1 << (m + 1), (1 << 0) | (1 << 2), (1 << 1) | (1 << 2)]
    matrices, physical_passes = [], 0
    for top0, top1 in itertools.product(range(4), repeat=2):
        if (((top0 >> 1) ^ top1) & 1) != ((nf['Gamma'][0] >> 1) & 1):
            continue
        ys = [nf['forward'][top0 | nf['offsets'][0]],
              nf['forward'][top1 | nf['offsets'][1]]]
        standard_columns = xs + ys
        columns = [lift_coordinate(coordinates(e, model.R + std) >> 1, standard_columns)
                   for e in model.E]
        assert model.criterion(columns) and model.direct_property(columns)
        physical_passes += 1
        matrices.append((top0, top1))
    assert physical_passes == 8
    diagonal_zero = [(u, v) for u, v in matrices if not (u & 1) and not (v >> 1 & 1)]
    assert len(diagonal_zero) == 2
    assert diagonal_zero[0][0] ^ diagonal_zero[1][0] == 2
    assert diagonal_zero[0][1] ^ diagonal_zero[1][1] == 1
    # Generic signed-edge control: an odd triangle has no orientations although
    # every single edge supplies a partner for each endpoint orientation.
    for signs in itertools.product(range(2), repeat=3):
        solutions = [v for v in itertools.product(range(2), repeat=3)
                     if (v[0] ^ v[1], v[1] ^ v[2], v[0] ^ v[2]) == signs]
        assert bool(solutions) == (sum(signs) % 2 == 0)
        assert all(any((u ^ v) == sign for v in range(2))
                   for sign in signs for u in range(2))
    return dict(status='PASS', quantum_parameters=[5, 1, 2], fixed_Z=xs,
                compatible_maps=8, direct_physical_passes=physical_passes,
                diagonal_zero_maps=2, phase_equivalence_classes=1,
                signed_triangle_controls=8,
                scope='Small real extension control and generic parity-cycle algebra; '
                      'not a new code construction or evidence of one-orbit domains at m9')


if __name__ == '__main__':
    report = verify()
    Path(__file__).with_name('one_orbit_control.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))
