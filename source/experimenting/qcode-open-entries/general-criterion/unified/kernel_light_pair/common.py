"""Experimental paired-light bound. Not a production producer."""
import json
import subprocess
import tempfile
from pathlib import Path
from rankone_residual import common as base
from rankone_light_budget.common import physical_data

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def forced(xs, entries, m, cut, purity):
    score = 0
    for a, weight in entries:
        zeros = sum(not ((a & x).bit_count() % 2) for x in xs)
        max_x = m - zeros
        max_x -= max_x % 2
        total = m + max_x + 2 * weight
        if total < 2 * purity:
            return None
        score += 2 if total < cut + purity else (1 if total < 2 * cut else 0)
    return score


def native(model, data, context, reps, *, backend, light=None, budget=None,
           group=True, fixed=None, max_nodes=None, seconds=None):
    command = [str((ROOT / ('kernel_z_compact/z_enum_compact' if backend == 'old-compact' else
                           'kernel_z_joint/z_enum_joint')) if backend.startswith('old-') else
                   HERE / f'z_{backend}'), '--gamma', '--avoid-radical', '--first-orbits', '--kernel-data']
    if backend.endswith('joint'):
        command += ['--mitm4', '65536']
    if fixed is not None:
        command += ['--prefix', ','.join(map(str, fixed))]
    if seconds is not None:
        command += ['--seconds', str(seconds)]
    if max_nodes is not None:
        command += ['--max-nodes', str(max_nodes)]
    text = base.payload(model, data, reps, context)
    if light is not None:
        assert budget is not None and budget >= 0
        command += ['--pair-budget', str(budget), '--pair-cut', str(light['cut']),
                    '--pair-purity', str(model.d)]
        text += str(len(light['entries'])) + '\n' + ''.join(f'{a} {w}\n' for a, w in light['entries'])
    with tempfile.TemporaryDirectory() as temporary:
        if group:
            path = Path(temporary) / 'group.txt'
            path.write_text(f'{len(context["group"])} {model.m - 1}\n' +
                            '\n'.join(' '.join(map(str, a)) for a in context['group']) + '\n')
            command += ['--prefix-group', str(path)]
        result = subprocess.run(command, input=text, text=True, capture_output=True, check=True, timeout=40)
    stats = json.loads(result.stderr)
    leaves = [tuple(map(int, line.split()[1:])) for line in result.stdout.splitlines()]
    assert stats['leaves'] == len(leaves) == len(set(leaves))
    return stats, leaves
