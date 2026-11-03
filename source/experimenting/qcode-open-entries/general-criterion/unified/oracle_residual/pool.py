"""Independent, disjoint Oracle shard; no reinterpretation of legacy hashes."""
import json
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'m8r1')]
from oracle_residual.common import HERE, assignment, fingerprints, activate, validate
from rankone_residual import common as base
from rankone_residual.search import expected_slices
from rankone_residual.aggregate import verify_chunk
from mixed_frames.stage3.common import atomic_json


def work(task):
    activate()
    row, s, k, hashes = task
    path = HERE/'runs'/s['names'][k]
    assert not path.exists(), path
    backend = 'joint' if s['prefix']['profile_index'] in (0, 8) else 'compact'
    r, ss = base.run_task(row, s['branch'], s['oi'], s['prefix'], s['ci'], k,
        None, hashes, candidate_dir=HERE/'candidates', backend=backend)
    r['oracle_producer'] = hashes['oracle_residual/verification.json']
    assert ss and r['orbit_digest'] == s['orbit_digest'] and r['chunks'] == s['chunks']
    atomic_json(path, r)
    return dict(task=path.name, status=r['status'], seconds=r['seconds'],
                backend=backend, nodes=r['z']['nodes'], leaves=r['z']['leaves'])


def main():
    assert __debug__
    plan = validate()
    hashes = fingerprints()
    assignment_data = assignment()
    rows = {r['cls']: r for r in plan['records']}
    (HERE/'runs').mkdir(exist_ok=True)
    for path in (HERE/'runs').glob('*.json'):
        assert json.loads(path.read_text())['status'] != 'FOUND', path

    def tasks():
        for cls in assignment_data['remote']:
            row = rows[cls]
            for s in expected_slices(plan, row):
                for k, name in enumerate(s['names']):
                    assert not (ROOT/'scheduler/STOP_REQUESTED').exists(), 'user stop'
                    assert shutil.disk_usage(HERE).free > 750_000_000, 'low VM disk space; preserve checkpoints and collect'
                    path = HERE/'runs'/name
                    if path.exists():
                        assert verify_chunk(path, row, s, k, hashes) == 'EXCLUDED_SLICE_CHUNK'
                        continue
                    yield row, s, k, hashes

    iterator = iter(tasks())
    with ProcessPoolExecutor(max_workers=2) as executor:
        pending = set()

        def submit():
            try:
                task = next(iterator)
            except StopIteration:
                return False
            pending.add(executor.submit(work, task))
            return True

        for _ in range(4):
            if not submit():
                break
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                r = future.result()
                print(json.dumps(r), flush=True)
                assert r['status'] == 'EXCLUDED_SLICE_CHUNK', 'FOUND/LIMIT remains unresolved; inspect'
                submit()
    atomic_json(HERE/'pool_finished.json', dict(status='BATCH_FINISHED_REQUIRES_STRICT_AGGREGATION', hashes=hashes))


if __name__ == '__main__':
    main()
