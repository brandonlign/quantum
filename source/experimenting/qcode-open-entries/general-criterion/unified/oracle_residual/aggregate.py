"""Strict union of legacy/Mac reuse and disjoint admitted Linux records."""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'m8r1')]
from oracle_residual.common import HERE, assignment, fingerprints
from rankone_reuse import aggregate as original
from rankone_reuse.records import verify_record as mac_verify, accepted_hashes
from rankone_residual.aggregate import verify_chunk
from mixed_frames.stage3.common import atomic_json, sha


def verify(path, row, s, k, *, accepted=None, require_pool=True):
    r = json.loads(Path(path).read_text())
    if r['hashes'] == accepted[2]:
        assert Path(path).parent.resolve() == (HERE/'runs').resolve()
        assert r.get('oracle_producer') == accepted[2]['oracle_residual/verification.json']
        assert not any(key in r for key in ('z_reuse', 'reuse_pool_hashes', 'pool_hashes'))
        return verify_chunk(Path(path), row, s, k, accepted[2])
    return mac_verify(path, row, s, k, accepted=accepted[:2], require_pool=require_pool)


def aggregate(plan, selected=None):
    a = assignment()
    selected = set(a['local']+a['remote']) if selected is None else set(selected)
    assert selected <= set(a['local']+a['remote'])
    accepted = accepted_hashes()+[fingerprints()]
    # This new process adapts the frozen aggregator's verifier, without changing
    # its source or the producer provenance in a single existing record.
    saved = original.verify_record
    original.verify_record = verify
    try:
        local = original.aggregate(plan, selected & set(a['local']), accepted=accepted)
        remote = original.aggregate(plan, selected & set(a['remote']),
            legacy=HERE/'empty_legacy', folder=HERE/'runs', accepted=accepted)
    finally:
        original.verify_record = saved
    errors = local['errors']+remote['errors']
    for directory, allowed in [(original.LEGACY, set(a['local'])),
                               (original.RUNS, set(a['local'])),
                               (HERE/'runs', set(a['remote']))]:
        for p in directory.glob('*.json'):
            try:
                assert int(p.name.split('_')[0]) in allowed, 'misassigned producer class'
            except Exception as e:
                errors.append(dict(path=str(p), error=repr(e)))
    records = sorted(local['records']+remote['records'], key=lambda r:r['cls'])
    orphans = local['orphans']+remote['orphans']
    complete = (len(selected) == len(plan['records']) and not errors and not orphans
        and len(records) == len(selected) and all(r['status'] ==
            'CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES' for r in records))
    return dict(status='CONDITIONALLY_COMPLETE_RESIDUAL_RANKONE' if complete else
            ('ERROR' if errors or orphans else 'OPEN'), classes=len(plan['records']),
        selected_classes=len(selected), counts=dict(Counter(local['counts'])+Counter(remote['counts'])),
        records=records, errors=errors, orphans=orphans,
        record_hashes=dict(local['record_hashes'], **remote['record_hashes']),
        accepted_hashes=accepted, aggregation_source_sha256=sha(Path(__file__)),
        assignment_sha256=sha(HERE/'assignment.json'),
        scope='Disjoint union of residual q1/q2 only; other branches and proof review remain separate.')


def main():
    from rankone_reuse.pool import validate
    ap = argparse.ArgumentParser()
    ap.add_argument('--classes', nargs='+', type=int)
    ap.add_argument('--output', type=Path, default=HERE/'aggregate.json')
    args = ap.parse_args()
    out = aggregate(validate(), args.classes)
    atomic_json(args.output, out)
    print(json.dumps({k:v for k,v in out.items() if k not in
        ('records','record_hashes','accepted_hashes')}), flush=True)


if __name__ == '__main__':
    main()
