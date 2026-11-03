"""Six-worker exhaustive first-image pipeline. Re-run to resume unfinished partitions.

The enumerator waits for OK after each leaf: at most one computation runs per
worker. A leaf is acknowledged only after exhaustive fixed_z completion and
the physical mu_out filter. Atomic checkpoints never promote LIMIT to COMPLETE.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
from pathlib import Path
import subprocess
import time
import traceback

from mixed_frames.fixed_z import solve
from .common import HERE, BINARY, atomic_json, graph_minimum, hashes, input_text, sha, target


def complete_leaf(model, frame, xs, mu):
    result = solve(model, frame, xs, find_all=True)
    assert result['status'] == 'COMPLETE' and result['complete'], result['status']
    accepted, minima = [], []
    for solution in result['solutions']:
        measured = graph_minimum(model, frame, solution['adapted_columns'])
        minima.append(measured)
        if measured >= mu:
            solution['minimum_stabilizer_outside_R'] = measured
            accepted.append(solution)
    return dict(X=xs, status='COMPLETE', complete=True, mu_out=mu,
                raw_solutions=len(result['solutions']), graph_minima=minima,
                accepted=accepted, completion=result)


def partition(task):
    row, first, expected_hashes = task
    folder = HERE/'runs'/f'{row["frame"]}_{row["class_index"]}'/str(first)
    folder.mkdir(parents=True, exist_ok=True)
    checkpoint = folder/'checkpoint.json'
    context = dict(frame=row['frame'], class_index=row['class_index'], first_image=first,
                   d=6, mu_out=row['mu_out'], input_sha256=row['input_sha256'], hashes=expected_hashes)
    if checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        assert all(saved[k] == v for k, v in context.items()), 'Checkpoint provenance changed'
        if saved['status'] == 'COMPLETE':
            return saved
    assert sha(HERE/row['input']) == row['input_sha256']
    model, frame = target(row['frame'], row['class_index'])
    completions = folder/'completions'
    completions.mkdir(exist_ok=True)
    count, rejected, surviving = 0, 0, 0
    start = time.monotonic()
    proc = None
    def save(status, **extra):
        value = dict(**context, status=status, complete=status == 'COMPLETE',
                     excluded=status == 'COMPLETE' and surviving == 0,
                     leaves=count, rejected_leaves=rejected, surviving_solutions=surviving,
                     seconds=time.monotonic()-start, pid=os.getpid(), **extra)
        atomic_json(checkpoint, value)
        return value
    save('RUNNING')
    try:
        with (folder/'stderr.log').open('w') as err, (folder/'enumerator.jsonl').open('w') as log:
            proc = subprocess.Popen([str(BINARY), str(HERE/row['input']), '--prefix', str(first), '--ack', '1'],
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err, text=True)
            summary = None
            for line in proc.stdout:
                log.write(line); log.flush()
                obj = json.loads(line)
                if obj['type'] == 'summary':
                    summary = obj
                    continue
                assert obj['type'] == 'leaf' and len(obj['X']) == frame.p
                xs = obj['X']
                count += 1
                key = sha_bytes(json.dumps(xs).encode())
                file = completions/(key+'.json')
                if file.exists():
                    done = json.loads(file.read_text())
                    assert done['X'] == xs and done['hashes'] == expected_hashes
                else:
                    done = complete_leaf(model, frame, xs, row['mu_out'])
                    done['hashes'] = expected_hashes
                    atomic_json(file, done)
                assert done['complete'] and done['status'] == 'COMPLETE'
                surviving += len(done['accepted'])
                rejected += not done['accepted']
                save('RUNNING', last_X=xs)
                proc.stdin.write('OK\n'); proc.stdin.flush()
            returncode = proc.wait()
            assert returncode == 0 and summary is not None, (returncode, summary)
            assert summary['leaf_count'] == count
            status = 'COMPLETE' if summary['complete'] else 'LIMIT'
            return save(status, enumerator=summary, log_sha256=sha(folder/'enumerator.jsonl'),
                        completion_sha256={p.name: sha(p) for p in sorted(completions.glob('*.json'))})
    except BaseException:
        if proc is not None and proc.poll() is None:
            proc.terminate(); proc.wait()
        return save('ERROR', error=traceback.format_exc())


def sha_bytes(data):
    import hashlib
    return hashlib.sha256(data).hexdigest()


def aggregate(verify=False):
    manifest = json.loads((HERE/'manifest.json').read_text())
    expected = manifest['hashes']
    assert hashes() == expected, 'Implementation/input source hashes changed'
    capacity = json.loads((HERE.parent/'capacity_certificate.json').read_text())['records']
    pairs = {(r['frame'],r['class_index']) for r in capacity if r['verdict']=='CANDIDATE'}
    assert len(manifest['classes']) == len(pairs) == 152
    assert {(r['frame'],r['class_index']) for r in manifest['classes']} == pairs
    classes = []
    for row in manifest['classes']:
        assert row['mu_out'] == (8 if row['frame'] in ('668','778') else 6)
        if verify:
            model, frame = target(row['frame'],row['class_index'])
            assert row['d'] == model.d == 6
            for k,v in dict(R=frame.R,E=frame.E,C=model.C,P=model.P,Gamma=model.Gamma,hull=model.hull).items():
                assert row[k] == v
            assert row['first_images'] == [s for s in range(1 << model.b)
                                           if frame.costs[1]+model.distances[s] >= model.d]
            assert (HERE/row['input']).read_text() == input_text(model,frame,row['mu_out'])
            assert sha(HERE/row['input']) == row['input_sha256']
        records = []
        for first in row['first_images']:
            path = HERE/'runs'/f'{row["frame"]}_{row["class_index"]}'/str(first)/'checkpoint.json'
            if not path.exists():
                records.append(dict(first_image=first, status='PENDING', complete=False, excluded=False))
                continue
            r = json.loads(path.read_text())
            assert r['hashes'] == expected and r['input_sha256'] == row['input_sha256']
            assert (r['frame'],r['class_index'],r['first_image']) == (row['frame'],row['class_index'],first)
            assert r['complete'] == (r['status'] == 'COMPLETE')
            assert r['excluded'] == (r['complete'] and r['surviving_solutions'] == 0)
            if r['complete']:
                assert r['enumerator']['complete'] and r['leaves'] == r['enumerator']['leaf_count']
                assert r['leaves'] == len(r['completion_sha256'])
                if verify:
                    assert sha(path.parent/'enumerator.jsonl') == r['log_sha256']
                    emitted = [json.loads(line)['X'] for line in (path.parent/'enumerator.jsonl').read_text().splitlines()
                               if json.loads(line)['type'] == 'leaf']
                    got = []
                    for name, digest in r['completion_sha256'].items():
                        cp = path.parent/'completions'/name
                        assert sha(cp) == digest
                        c = json.loads(cp.read_text())
                        assert c['complete'] and c['status'] == 'COMPLETE' and c['hashes'] == expected
                        assert c['completion']['complete'] and c['completion']['status'] == 'COMPLETE'
                        assert c['X'] == c['completion']['X'] and c['mu_out'] == row['mu_out']
                        assert c['raw_solutions'] == len(c['completion']['solutions'])
                        assert c['graph_minima'] == [graph_minimum(model,frame,s['adapted_columns'])
                                                     for s in c['completion']['solutions']]
                        assert c['accepted'] == [s for s,minimum in zip(c['completion']['solutions'],c['graph_minima'])
                                                if minimum >= row['mu_out']]
                        got.append(c['X'])
                    assert sorted(emitted) == sorted(got) and len({tuple(x) for x in got}) == len(got)
            r['checkpoint_sha256'] = sha(path)
            records.append(r)
        complete = all(r['complete'] for r in records)
        status = 'COMPLETE' if complete else ('ERROR' if any(r['status']=='ERROR' for r in records) else 'INCOMPLETE')
        classes.append(dict(frame=row['frame'], class_index=row['class_index'], mu_out=row['mu_out'],
                            hull=row['hull'], status=status, complete=complete,
                            excluded=complete and all(r['excluded'] for r in records),
                            partition_count=len(records), complete_partitions=sum(r['complete'] for r in records),
                            leaves=sum(r.get('leaves',0) for r in records), partitions=records))
    value = dict(hashes=expected, manifest_sha256=sha(HERE/'manifest.json'), classes=classes,
                 status='COMPLETE' if all(c['complete'] for c in classes) else 'INCOMPLETE',
                 verified=verify, updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                 scope='Conditional on supplied census and profile/frame certificates. LIMIT/ERROR/PENDING/RUNNING never excludes.')
    atomic_json(HERE/'aggregate.json', value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--aggregate', action='store_true')
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--frames', nargs='+', default=['668', '778', '688'])
    args = parser.parse_args()
    assert 1 <= args.workers <= 6
    if args.aggregate:
        r = aggregate(args.verify)
        print(json.dumps(dict(status=r['status'], complete_classes=sum(c['complete'] for c in r['classes']),
                              excluded_classes=sum(c['excluded'] for c in r['classes']))), flush=True)
        return
    manifest = json.loads((HERE/'manifest.json').read_text())
    assert hashes() == manifest['hashes']
    # Finish both mu=8 frames before starting 688.
    for frame in args.frames:
        rows = [r for r in manifest['classes'] if r['frame'] == frame]
        tasks = [(r, first, manifest['hashes']) for r in rows for first in r['first_images']]
        print(json.dumps(dict(frame=frame, partitions=len(tasks), workers=args.workers)), flush=True)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(partition, t) for t in tasks]
            for future in as_completed(futures):
                r = future.result()
                print(json.dumps({k:r[k] for k in ('frame','class_index','first_image','status','leaves','seconds')}), flush=True)
                aggregate()
    aggregate(True)


if __name__ == '__main__':
    main()
