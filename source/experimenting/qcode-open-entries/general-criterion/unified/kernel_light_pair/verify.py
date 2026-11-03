"""Independent physical pair bounds and complete physical small-model maps."""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from kernel_light_pair.common import forced, native, physical_data
from model import ExtensionModel, dual, lift_coordinate, span, split, weight
from gauge_fix import supported_subspace
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts
from fixed_kernel.solve import solve
from mixed_frames.stage3.common import atomic_json, sha


def main():
    started = time.monotonic()
    physical_pairs = 0
    # The bound is computed from X bits; the oracle uses both physical words.
    for m, d, cut in [(4, 3, 4), (6, 3, 5), (8, 6, 8)]:
        R = ((1 << m) - 1) << m
        for h in range(1 << (2 * m)):
            x = h & ((1 << m) - 1)
            if x.bit_count() % 2:
                continue
            for suffix_weight in range(1, 4):
                pair = [weight(h, m) + suffix_weight, weight(h ^ R, m) + suffix_weight]
                if min(pair) < d:
                    continue
                actual = sum(w < cut for w in pair)
                for depth in range(m):
                    bits = [(x >> i) & 1 for i in range(depth)]
                    bound = forced(bits, [(1, suffix_weight)], m, cut, d)
                    assert bound is not None and bound <= actual, (m, h, suffix_weight, depth, bound, pair)
                physical_pairs += 1
    fixtures = []
    a = json.loads((ROOT / 'fixed_kernel/control.json').read_text())
    fixtures.append((a['H'], ExtensionModel(8, 2, 3, 4, [240], a['C']), a['N']))
    S = json.loads((ROOT / 'm8r1/fixed_radical_verification.json').read_text())['nonzero_radical_positive_control']['generators']
    H = [split(w, 4, 3)[0] for w in supported_subspace(dual(S, 7), 7, range(4))]
    fixtures.append((H, ExtensionModel(7, 1, 3, 4, [240], [7]), []))
    S = json.loads((ROOT / 'radical_projection/positive_8.json').read_text())['result']['generators']
    D = dual(S, 8)
    H = [split(w, 4, 4)[0] for w in supported_subspace(D, 8, range(4))]
    C = [split(w, 4, 4)[1] for w in supported_subspace(D, 8, range(4, 8))]
    N = [split(w, 4, 4)[1] for w in supported_subspace(S, 8, range(4, 8))]
    fixtures.append((H, ExtensionModel(8, 2, 3, 4, [240], C), N))
    records = []
    positive_maps = 0
    for H, model, N in fixtures:
        data = quotient_input(model, N)
        H, order, contexts = kernel_contexts(H, model.m, data['Gamma'], with_group=True)
        lights = {cut: physical_data(model, N, data, cut) for cut in (4, 5, 6)}
        reps = [i for i, d in enumerate(data['distances']) if d >= model.d - 1]
        for ci, context in enumerate(contexts):
            old, leaves = native(model, data, context, reps, backend='old-compact', group=False)
            assert old['status'] == 'COMPLETE'
            scores = {cut: {xs: forced(xs, light['entries'], model.m, cut, model.d) for xs in leaves}
                      for cut, light in lights.items()}
            wide_controls = {}
            for backend in ['compact', 'joint']:
                unfiltered, unchanged = native(model, data, context, reps, backend=backend, group=False)
                assert unfiltered['status'] == 'COMPLETE' and unchanged == leaves
                for cut, light in lights.items():
                    for budget in [0, 1, 2, 10, 64]:
                        stats, filtered = native(model, data, context, reps, backend=backend, group=False,
                                                 light=light, budget=budget)
                        expected = [xs for xs in leaves if scores[cut][xs] is not None and scores[cut][xs] <= budget]
                        assert stats['status'] == 'COMPLETE' and filtered == expected
                        if budget == 64:
                            wide_controls[backend, cut] = set(filtered)
                        records.append(dict(parameters=[model.n, model.k, model.d], context=ci,
                                            backend=backend, cut=cut, budget=budget, leaves=len(leaves),
                                            retained=len(filtered), pruned=stats.get('pair_pruned', 0)))
            for xs in leaves:
                answer = solve(model, H, [lift_coordinate(x, data['units']) for x in xs], find_all=True)
                assert answer['status'] == 'COMPLETE'
                for solution in answer['solutions']:
                    S, D = model.construct(solution['model_columns'])
                    assert min(weight(w, model.n) for w in span(D) if w) >= model.d
                    weights = [weight(w, model.n) for w in span(S) if w]
                    for cut, light in lights.items():
                        actual = sum(w < cut for w in weights)
                        assert scores[cut][xs] is not None and scores[cut][xs] <= actual - light['baseline']
                        assert all(xs in wide_controls[backend, cut] for backend in ('compact', 'joint'))
                    positive_maps += 1
            for backend in ['compact', 'joint']:
                limited, _ = native(model, data, context, reps, backend=backend, group=False,
                                    light=lights[6], budget=10, max_nodes=0)
                assert limited['status'] == 'LIMIT'
    assert positive_maps > 0
    active = sum(r['pruned'] > 0 for r in records)
    assert active > 0, 'Native active-rejection path must be exercised'
    files = [Path(__file__), HERE / 'common.py', HERE / 'PROOF.md', HERE / 'build.py', HERE / 'build.json',
             ROOT / 'rankone_light_budget/common.py', ROOT / 'model.py', ROOT / 'fixed_kernel/solve.py']
    files += [HERE / f'z_{backend}{suffix}' for backend in ('compact', 'joint') for suffix in ('', '.c')]
    report = dict(status='PASS_BOUNDED_PAIR_PRUNE_CONTROLS', physical_pair_partial_bounds=physical_pairs,
                  genuine_positive_maps_checked=positive_maps, complete_comparisons=records,
                  active_rejection_comparisons=active,
                  hashes={str(p.relative_to(ROOT)): sha(p) for p in files}, seconds=time.monotonic() - started,
                  scope='Experimental prune correctness controls only; no target exclusion or production deployment')
    atomic_json(HERE / 'verification.json', report)
    print(json.dumps({k: v for k, v in report.items() if k not in ('hashes', 'complete_comparisons')}), flush=True)


if __name__ == '__main__':
    main()
