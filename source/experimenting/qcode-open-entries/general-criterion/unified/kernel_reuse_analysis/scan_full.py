"""Permit dependent-site motion; still fixes the first label. Analysis only."""
import itertools
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from mixed_frames.stage3.common import atomic_json, sha


def lift(mask, rows):
    out = 0
    for i, row in enumerate(rows):
        if mask & (1 << i):
            out ^= row
    return out


def matrices(order):
    inv = [order.index(i) for i in range(8)]
    return ([1 << order[i] if order[i] < 7 else 127 for i in range(7)],
            [1 << inv[i] if inv[i] < 7 else 127 for i in range(7)])


def conjugate(action, P, inverse):
    return tuple(lift(lift(v, action), P) for v in inverse)


def congruence(gram, inverse):
    return tuple(sum(((lift(v, gram) & w).bit_count() % 2) << j
                     for j, w in enumerate(inverse)) for v in inverse)


def permutations(left, right):
    if left[0] != right[0] or sorted(left) != sorted(right):
        return
    colors = sorted(set(left))
    source = {c: [i for i in range(1, 8) if left[i] == c] for c in colors}
    target = {c: [i for i in range(1, 8) if right[i] == c] for c in colors}
    for choices in itertools.product(*(itertools.permutations(target[c]) for c in colors)):
        order = list(range(8))
        for c, images in zip(colors, choices):
            for i, j in zip(source[c], images):
                order[i] = j
        yield order


def main():
    original = json.loads((HERE / 'scan.json').read_text())
    results = []
    for case in original['cases']:
        contexts = case['contexts']
        arrows = []
        for i, a in enumerate(contexts):
            for j, b in enumerate(contexts):
                if i == j or a['group_order'] > b['group_order']:
                    continue
                target = set(map(tuple, b['group']))
                for order in permutations(a['tags'], b['tags']):
                    P, inverse = matrices(order)
                    assert all(lift(P[k], inverse) == 1 << k for k in range(7))
                    if congruence(a['Gram'], inverse) != tuple(b['Gram']):
                        continue
                    if all(conjugate(A, P, inverse) in target for A in a['group']):
                        arrows.append(dict(source=i, target=j, source_profile=a['profile'], target_profile=b['profile'],
                                           source_identification=a['identification'], target_identification=b['identification'],
                                           site_image=order, source_group_order=a['group_order'], target_group_order=b['group_order']))
                        break
        result = dict(q=case['q'], t=case['t'], contexts=len(contexts), containments=arrows)
        results.append(result)
        print(json.dumps(result), flush=True)
    atomic_json(HERE / 'scan_full.json', dict(status='BOUNDED_STRUCTURAL_ANALYSIS_ONLY', cases=results,
                hashes={str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), HERE / 'scan.json']},
                scope='Potential input conjugacies only; no coverage claim, deployment or target exclusion.'))


if __name__ == '__main__':
    main()
