"""Resume legacy coverage, reusing only proven identical profile-1 Z streams."""
import argparse
import json
import sys
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from rankone_residual import common as baseline
from rankone_residual.search import expected_slices
from rankone_reuse.records import LEGACY, RUNS, accepted_hashes, pool_hashes, record_path, verify_record
from mixed_frames.stage3.common import atomic_json, sha


def validate():
    ready = json.loads((HERE / 'READY.json').read_text())
    assert ready['status'] == 'READY_FOR_RESIDUAL_REUSE_SWEEP'
    for path, digest in ready['hashes'].items():
        assert sha(ROOT / path) == digest, path
    return baseline.validate(require_controls=True)


def work(task):
    from rankone_reuse.common import run_task
    row, s, k, hashes, provenance, folder, legacy, backend = task
    folder, legacy = Path(folder), Path(legacy)
    path = folder / s['names'][k]
    assert not path.exists(), ('overlapping dispatch', path)
    source = None
    if (s['q'], s['t'], s['prefix']['profile_index'], s['ci'], backend) == (1, 0, 2, 0, 'compact'):
        name = baseline.task_name(row['cls'], 1, 0, s['oi'], 1, 0, k)
        source = record_path(name, legacy=legacy, folder=folder)
        if source is not None:
            status = json.loads(source.read_text())['status']
            if status == 'FOUND':
                raise RuntimeError(f'Inspect target candidate {source}')
            if status != 'EXCLUDED_SLICE_CHUNK':
                source = None
    r = run_task(row, s, k, hashes, source_path=source,
                 candidate_dir=folder / 'candidates', backend=backend)
    if provenance is not None:
        r['reuse_pool_hashes'] = provenance
    assert r['chunks'] == s['chunks'] and r['orbit_digest'] == s['orbit_digest']
    atomic_json(path, r)
    return dict(task=path.name, status=r['status'], backend=backend,
                reused='z_reuse' in r, seconds=r['seconds'],
                nodes=r['z']['nodes'], leaves=r['z']['leaves'])


def work_unit(tasks):
    """Finish a new source before its paired target, without waiting workers."""
    results = []
    for task in tasks:
        result = work(task)
        results.append(result)
        if result['status'] != 'EXCLUDED_SLICE_CHUNK':
            break
    return results


def main():
    from rankone_reuse.common import fingerprints
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--classes', type=int, nargs='+')
    ap.add_argument('--profiles', type=int, nargs='+')
    ap.add_argument('--max-tasks', type=int)
    ap.add_argument('--backend', choices=['joint-auto'], default='joint-auto')
    args = ap.parse_args()
    assert __debug__ and 1 <= args.workers <= 8
    assert args.max_tasks is None or args.max_tasks > 0
    plan = validate()
    hashes, accepted, provenance = fingerprints(), accepted_hashes(), pool_hashes()
    RUNS.mkdir(exist_ok=True)
    selected = None if args.classes is None else set(args.classes)
    assert selected is None or selected <= {r['cls'] for r in plan['records']}
    # Fail closed even for a candidate outside the requested test selection.
    for directory in (LEGACY, RUNS):
        for path in directory.glob('*.json'):
            if json.loads(path.read_text()).get('status') == 'FOUND':
                raise RuntimeError(f'Inspect saved target candidate {path}')

    def tasks():
        count = 0
        for row in plan['records']:
            if selected is not None and row['cls'] not in selected:
                continue
            scheduled = set()
            companions = {(s['oi'], s['ci']): s for s in
                          expected_slices(plan, row, profiles=[2])}
            for s in expected_slices(plan, row, profiles=args.profiles):
                for k, name in enumerate(s['names']):
                    if name in scheduled:
                        continue
                    assert not (ROOT / 'scheduler/STOP_REQUESTED').exists(), 'user stop flag'
                    path = record_path(name)
                    if path is not None:
                        status = verify_record(path, row, s, k, accepted=accepted)
                        if status == 'FOUND':
                            raise RuntimeError(f'Inspect target candidate {path}')
                        if status == 'EXCLUDED_SLICE_CHUNK':
                            continue
                        # Preserve unresolved versioned evidence rather than overwrite it.
                        assert path.parent != RUNS, ('unresolved reuse record', path)
                    backend = 'joint' if s['prefix']['profile_index'] in (0, 8) else 'compact'
                    unit = [(row, s, k, hashes, provenance, str(RUNS), str(LEGACY), backend)]
                    # Bind newly computed profile-1 and profile-2 chunks to one
                    # worker so eight-way execution cannot race past its source.
                    if (s['q'], s['t'], s['prefix']['profile_index'], s['ci']) == (1, 0, 1, 0) and \
                       (args.profiles is None or 2 in args.profiles) and \
                       (args.max_tasks is None or count + 2 <= args.max_tasks):
                        partner = companions[s['oi'], 0]
                        partner_name = partner['names'][k]
                        assert partner['chunks'] == s['chunks']
                        old = record_path(partner_name)
                        status = None if old is None else verify_record(old, row, partner, k, accepted=accepted)
                        if status == 'FOUND':
                            raise RuntimeError(f'Inspect target candidate {old}')
                        if status != 'EXCLUDED_SLICE_CHUNK':
                            assert old is None or old.parent != RUNS, ('unresolved reuse record', old)
                            unit.append((row, partner, k, hashes, provenance, str(RUNS), str(LEGACY), 'compact'))
                            scheduled.add(partner_name)
                    yield unit
                    count += len(unit)
                    if args.max_tasks is not None and count >= args.max_tasks:
                        return

    iterator = iter(tasks())
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        pending = set()

        def submit():
            try:
                task = next(iterator)
            except StopIteration:
                return False
            pending.add(executor.submit(work_unit, task))
            return True

        for _ in range(2 * args.workers):
            if not submit():
                break
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                for result in future.result():
                    print(json.dumps(result), flush=True)
                    if result['status'] == 'FOUND':
                        for f in pending:
                            f.cancel()
                        raise RuntimeError('Candidate persisted; stop dispatch and inspect it')
                    assert result['status'] == 'EXCLUDED_SLICE_CHUNK', 'Unfinished chunk stays unresolved'
                submit()
    atomic_json(HERE / 'pool_finished.json', dict(status='BATCH_FINISHED_REQUIRES_STRICT_AGGREGATION',
                args=vars(args), hashes=hashes, reuse_pool_hashes=provenance))


if __name__ == '__main__':
    main()
