"""Bounded native/reuse pool equivalence and strict mixed-producer checks."""
import copy
import json
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from rankone_residual import common as base
from rankone_residual.search import expected_slices
from rankone_reuse import common
from rankone_reuse.aggregate import aggregate
from rankone_reuse.pool import work, work_unit
from rankone_reuse.records import LEGACY, verify_record, record_path
from mixed_frames.stage3.common import atomic_json, sha


def scrub(x):
    if isinstance(x, dict):
        return {k: scrub(v) for k, v in x.items() if k != 'seconds'}
    if isinstance(x, list):
        return [scrub(v) for v in x]
    return x


def main():
    plan = base.validate(require_controls=True)
    hashes = common.fingerprints()
    accepted = [base.fingerprints(), hashes]
    lookup = {r['cls']: r for r in plan['records']}
    row = lookup[0]
    slices = list(expected_slices(plan, row))
    by_profile = {s['prefix']['profile_index']: s for s in slices}
    assert {0, 1, 2, 3, 8, 12} <= set(by_profile)
    control = HERE / 'integration_control'
    assert not control.exists(), 'Preserve previous verification; use a new version for changes'
    control.mkdir()
    legacy = control / 'legacy'
    legacy.mkdir()
    for s in slices:
        for name in s['names']:
            shutil.copyfile(LEGACY / name, legacy / name)
    baseline_report = aggregate(plan, [0], legacy=legacy, folder=control / 'empty',
                                accepted=accepted, require_pool=False)
    assert baseline_report['status'] == 'OPEN' and baseline_report['counts']['CLOSED_CLASSES'] == 1
    chosen = [(row, by_profile[2], 0), (row, by_profile[2], 1),
              (row, by_profile[0], 0), (row, by_profile[3], 0), (row, by_profile[8], 0)]
    row5 = lookup[5]
    s5 = next(s for s in expected_slices(plan, row5, profiles=[4]) if s['oi'] == 0)
    chosen.append((row5, s5, 0))
    for r, s, k in chosen:
        if r['cls'] == 0:
            (legacy / s['names'][k]).unlink()
    folders = [control / 'sequential', control / 'parallel']
    for folder in folders:
        folder.mkdir()
    tasks = [[(r, s, k, hashes, None, str(folder), str(legacy),
               'joint' if s['prefix']['profile_index'] in (0, 8) else 'compact')
              for r, s, k in chosen] for folder in folders]
    a = [work(t) for t in tasks[0]]
    with ProcessPoolExecutor(max_workers=2) as executor:
        b = list(executor.map(work, tasks[1]))
    assert all(x['status'] == 'EXCLUDED_SLICE_CHUNK' for x in a + b)
    assert sum(x['reused'] for x in a) == sum(x['reused'] for x in b) == 2
    for r, s, k in chosen:
        name = s['names'][k]
        x, y = [json.loads((f / name).read_text()) for f in folders]
        assert scrub(x) == scrub(y), name
        verify_record(folders[0] / name, r, s, k, accepted=accepted, require_pool=False)
        original = json.loads((LEGACY / name).read_text())
        # Baseline solver/backend and each H's own Y evidence must agree.
        for field in ('z', 'leaves', 'tally', 'H', 'first', 'status', 'backend'):
            assert scrub(x[field]) == scrub(original[field]), (name, field)
    reports = [aggregate(plan, [0], legacy=legacy, folder=f,
                         accepted=accepted, require_pool=False) for f in folders]
    for report in reports:
        assert report['status'] == 'OPEN' and not report['errors'] and not report['orphans']
        assert report['counts']['CLOSED_CLASSES'] == 1
        assert report['counts']['REUSED_Z_CHUNKS'] == 2
        assert {k: v for k, v in report['counts'].items() if k != 'REUSED_Z_CHUNKS'} == baseline_report['counts']
    path = folders[0] / by_profile[2]['names'][0]
    original = json.loads(path.read_text())
    rejected = []

    def reject(label, value, *, require_pool=False):
        atomic_json(path, value)
        try:
            verify_record(path, row, by_profile[2], 0, accepted=accepted, require_pool=require_pool)
        except (AssertionError, KeyError, ValueError, RuntimeError, FileNotFoundError):
            rejected.append(label)
        else:
            raise AssertionError(('corruption accepted', label))
        finally:
            atomic_json(path, original)

    reject('missing production pool provenance', original, require_pool=True)
    bad = copy.deepcopy(original); bad['hashes'] = {'unknown': 'producer'}; reject('unknown producer', bad)
    bad = copy.deepcopy(original); bad['hashes'] = accepted[0]; reject('legacy producer hiding reuse', bad)
    bad = copy.deepcopy(original); bad['z']['status'] = 'LIMIT'; reject('Z LIMIT relabeled as exclusion', bad)
    bad = copy.deepcopy(original); bad['leaves'] = []; reject('missing Y leaves', bad)
    bad = copy.deepcopy(original); bad['z_reuse']['source_sha256'] = '0' * 64; reject('changed source digest', bad)
    source = ROOT / original['z_reuse']['source_path']
    saved_source = source.read_bytes()
    try:
        bad_source = json.loads(saved_source); bad_source['status'] = 'LIMIT'
        atomic_json(source, bad_source)
        reject('source became unresolved', original)
    finally:
        source.write_bytes(saved_source)
    bad = copy.deepcopy(original); bad['status'] = 'FOUND'; atomic_json(path, bad)
    found_report = aggregate(plan, [5], legacy=legacy, folder=folders[0], accepted=accepted, require_pool=False)
    assert found_report['status'] == 'ERROR' and found_report['errors']
    rejected.append('FOUND outside selected class blocks completion')
    atomic_json(path, original)
    with tempfile.TemporaryDirectory(dir=control) as temp:
        duplicate = Path(temp)
        name = by_profile[1]['names'][0]
        shutil.copyfile(legacy / name, duplicate / name)
        try:
            record_path(name, legacy=legacy, folder=duplicate)
        except AssertionError:
            rejected.append('conflicting completed producers')
        else:
            raise AssertionError('conflicting producers accepted')
    paired = control / 'paired_unit'
    paired.mkdir()
    empty = control / 'no_legacy'
    empty.mkdir()
    pair_tasks = [(row, by_profile[p], 0, hashes, None, str(paired), str(empty), 'compact') for p in (1, 2)]
    with ProcessPoolExecutor(max_workers=2) as executor:
        pair_results = executor.submit(work_unit, pair_tasks).result()
    assert len(pair_results) == 2 and not pair_results[0]['reused'] and pair_results[1]['reused']
    for p in (1, 2):
        path2 = paired / by_profile[p]['names'][0]
        verify_record(path2, row, by_profile[p], 0, accepted=accepted, require_pool=False)
        record = json.loads(path2.read_text())
        original2 = json.loads((LEGACY / path2.name).read_text())
        for field in ('z', 'leaves', 'tally', 'status'):
            assert scrub(record[field]) == scrub(original2[field])
    report = dict(status='PASS', complete_chunks=len(chosen) + 2, paired_processes=2,
                  reused_per_run=2, target_native_evidence_matches=True,
                  new_producer_source_and_paired_target_verified=True,
                  strict_closed_class_unchanged=True, rejected=rejected,
                  hashes={str(p.relative_to(ROOT)): sha(p) for p in
                          [HERE / name for name in ('pool.py', 'records.py', 'aggregate.py', 'verify_integration.py')]},
                  engine_hashes=hashes)
    atomic_json(HERE / 'integration_verification.json', report)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
