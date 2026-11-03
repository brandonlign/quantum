"""Bounded target probes only; do not treat these as whole-case exclusions."""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from kernel_light_pair.common import forced, native, physical_data
from rankone_residual import common as base
from rankone_residual.search import expected_slices
from mixed_frames.stage3.common import atomic_json, sha


def main():
    proof = json.loads((HERE / 'verification.json').read_text())
    assert proof['status'] == 'PASS_BOUNDED_PAIR_PRUNE_CONTROLS'
    for name, digest in proof['hashes'].items():
        assert sha(ROOT / name) == digest, name
    plan = base.validate(require_controls=True)
    rows = {r['cls']: r for r in plan['records']}
    chosen = [(0, p) for p in (0, 1, 3, 8, 12)] + [(1, p) for p in (0, 8)] + [(5, p) for p in (4, 5, 6, 9, 13, 16)]
    results = []
    for cls, profile in chosen:
        row = rows[cls]
        s = next(s for s in expected_slices(plan, row, profiles=[profile]) if s['oi'] == 0)
        model, data, H, order, contexts, io = base.slice_setup(tuple(row['C']), tuple(s['branch']['N_orbits'][0]['N']),
                    tuple(s['prefix']['H']), s['prefix']['p'], s['q'], s['t'])
        c = contexts[s['ci']]
        ss = base.suffix_slice(model, s['branch']['N_orbits'][0]['N'], data, c)
        light = physical_data(model, s['branch']['N_orbits'][0]['N'], data, 8)
        assert light['baseline'] <= 2
        saved = json.loads((ROOT / 'rankone_residual/runs' / s['names'][0]).read_text())
        reps = sorted(set(ss['first_representatives'][:1] + [saved['leaves'][0]['X'][0]] if saved['leaves'] else ss['first_representatives'][:1]))
        backend = 'joint' if profile in (0, 8) else 'compact'
        measured = []
        for mode in ('old', 'pair'):
            start = time.monotonic()
            stats, leaves = native(model, data, c, reps, backend=('old-' if mode == 'old' else '') + backend,
                                   light=light if mode == 'pair' else None,
                                   budget=2 - light['baseline'] if mode == 'pair' else None, seconds=10)
            measured.append(dict(mode=mode, wall=time.monotonic() - start, stats=stats, leaves=leaves))
        if all(m['stats']['status'] == 'COMPLETE' for m in measured):
            expected = [xs for xs in measured[0]['leaves'] if forced(xs, light['entries'], 8, 8, 6) is not None
                        and forced(xs, light['entries'], 8, 8, 6) <= 2 - light['baseline']]
            assert measured[1]['leaves'] == expected
        result = dict(cls=cls, profile=profile, context=s['ci'], reps=reps, baseline=light['baseline'],
                      physical_entries=len(light['entries']), measurements=measured)
        results.append(result)
        print(json.dumps({k: v for k, v in result.items() if k != 'measurements'} | dict(
              measurements=[{k: v for k, v in m.items() if k != 'leaves'} for m in measured])), flush=True)
        atomic_json(HERE / 'benchmark.json', dict(status='BOUNDED_TARGET_PROBES_ONLY', results=results,
                    hashes={str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), HERE / 'verification.json']},
                    scope='Necessary-prune timing experiment only; no production deployment or target exclusion'))


if __name__ == '__main__':
    main()
