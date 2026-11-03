"""Bounded search for coordinate-conjugate native problems. Not production."""
import itertools
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from rankone_residual.common import validate, slice_setup
from rankone_residual.search import expected_slices
from mixed_frames.stage3.common import atomic_json, sha


def permutations(left, right):
    if left[0] != right[0] or left[7] != right[7]:
        return
    colors = sorted(set(left))
    source = {c: [i for i in range(1, 7) if left[i] == c] for c in colors}
    target = {c: [i for i in range(1, 7) if right[i] == c] for c in colors}
    if any(len(source[c]) != len(target[c]) for c in colors):
        return
    if sorted(left) != sorted(right):
        return
    for choices in itertools.product(*(itertools.permutations(target[c]) for c in colors)):
        order = list(range(8))
        for c, images in zip(colors, choices):
            for i, j in zip(source[c], images):
                order[i] = j
        assert sorted(order) == list(range(8))
        yield order


def transform(matrix, order):
    out = [0] * 7
    for i, row in enumerate(matrix):
        out[order[i]] = sum(((row >> j) & 1) << order[j] for j in range(7))
    return tuple(out)


def scan(plan, row, q, t):
    slices = [s for s in expected_slices(plan, row, branches={(q, t)}) if s['oi'] == 0]
    contexts = []
    for s in slices:
        prefix = s['prefix']
        model, data, H, order, all_contexts, io = slice_setup(
            tuple(row['C']), tuple(s['branch']['N_orbits'][0]['N']), tuple(prefix['H']), prefix['p'], q, t)
        c = all_contexts[s['ci']]
        contexts.append(dict(profile=prefix['profile_index'], identification=s['ci'], tags=c['tags'],
                             Gram=c['Gram'], group=c['group'], group_order=len(c['group'])))
    arrows = []
    for i, a in enumerate(contexts):
        for j, b in enumerate(contexts):
            if i == j:
                continue
            target = set(map(tuple, b['group']))
            if len(a['group']) > len(target):
                continue
            checked = 0
            for order in permutations(a['tags'], b['tags']):
                checked += 1
                if transform(a['Gram'], order) != tuple(b['Gram']):
                    continue
                if all(transform(action, order) in target for action in a['group']):
                    arrows.append(dict(source=i, target=j, source_profile=a['profile'], target_profile=b['profile'],
                                       source_identification=a['identification'], target_identification=b['identification'],
                                       site_image=order, source_group_order=len(a['group']), target_group_order=len(target),
                                       permutations_checked=checked))
                    break
    return dict(q=q, t=t, cls=row['cls'], contexts=contexts, proven_input_containments=arrows)


def main():
    start = time.monotonic()
    plan = validate(require_controls=True)
    rows = {r['cls']: r for r in plan['records']}
    reports = [scan(plan, rows[0], 1, 0), scan(plan, rows[0], 2, 2), scan(plan, rows[5], 2, 0)]
    for r in reports:
        print(json.dumps(dict(q=r['q'], t=r['t'], contexts=len(r['contexts']), containments=r['proven_input_containments'])), flush=True)
    atomic_json(HERE / 'scan.json', dict(status='BOUNDED_STRUCTURAL_ANALYSIS_ONLY', cases=reports,
                seconds=time.monotonic() - start,
                hashes={str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), ROOT / 'rankone_residual/READY.json',
                        ROOT / 'kernel_case_plan/plan.json', ROOT / 'fixed_kernel/common.py']},
                scope='No deployment, no native run, and no target exclusion. Potential full-slice reuse needs its own proof and independent validation.'))


if __name__ == '__main__':
    main()
