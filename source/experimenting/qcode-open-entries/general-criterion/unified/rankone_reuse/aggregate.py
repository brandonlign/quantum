"""Strict union of untouched legacy evidence and separately versioned reuse."""
import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from rankone_residual.common import validate
from rankone_residual.search import expected_slices
from rankone_reuse.records import LEGACY, RUNS, accepted_hashes, record_path, verify_record
from mixed_frames.stage3.common import atomic_json, sha


def aggregate(plan, selected, *, legacy=LEGACY, folder=RUNS, require_pool=True, accepted=None):
    lookup = {r['cls']: r for r in plan['records']}
    selected = set(lookup) if selected is None else set(selected)
    assert selected <= set(lookup)
    accepted = accepted_hashes() if accepted is None else accepted
    present = defaultdict(set)
    orphans, errors, records = [], [], []
    counts, record_hashes = Counter(), {}
    for directory in (Path(legacy), Path(folder)):
        for path in directory.glob('*.json'):
            match = re.fullmatch(r'(\d+)_q[12]t[02]_\d+_\d+_\d+_\d+\.json', path.name)
            if not match or int(match[1]) not in lookup:
                orphans.append(str(path))
                continue
            cls = int(match[1])
            # A positive result anywhere keeps the whole census unresolved.
            if json.loads(path.read_text()).get('status') == 'FOUND':
                errors.append(dict(path=str(path), error='Saved FOUND candidate requires inspection'))
            if cls in selected:
                present[cls].add(path.name)
    for cls in sorted(selected):
        row = lookup[cls]
        if not present[cls]:
            counts['PENDING_CLASSES'] += 1
            records.append(dict(cls=cls, status='PENDING', branches=[dict(q=b['q'], t=b['t'], status='PENDING') for b in row['branches']]))
            continue
        seen, branch_counts = set(), defaultdict(Counter)
        for s in expected_slices(plan, row):
            ct = Counter()
            for k, name in enumerate(s['names']):
                seen.add(name)
                try:
                    path = record_path(name, legacy=legacy, folder=folder)
                    if path is None:
                        ct['PENDING'] += 1
                        continue
                    status = verify_record(path, row, s, k, accepted=accepted, require_pool=require_pool)
                    ct[status] += 1
                    record_hashes[str(path.resolve())] = sha(path)
                    r = json.loads(path.read_text())
                    if 'z_reuse' in r:
                        counts['REUSED_Z_CHUNKS'] += 1
                        source = (ROOT / r['z_reuse']['source_path']).resolve()
                        record_hashes[str(source)] = sha(source)
                except Exception as e:
                    ct['ERROR'] += 1
                    errors.append(dict(cls=cls, task=name, error=repr(e)))
            counts.update(ct)
            closed = ct == Counter(EXCLUDED_SLICE_CHUNK=s['chunks'])
            branch_counts[s['q'], s['t']]['CLOSED_SLICES' if closed else 'OPEN_SLICES'] += 1
        orphans.extend(str(Path(folder) / name) for name in sorted(present[cls] - seen))
        branches = [dict(q=b['q'], t=b['t'],
                    status='CONDITIONALLY_EXCLUDED_BRANCH' if branch_counts[b['q'], b['t']]['CLOSED_SLICES'] > 0 and branch_counts[b['q'], b['t']]['OPEN_SLICES'] == 0 else 'OPEN',
                    counts=dict(branch_counts[b['q'], b['t']])) for b in row['branches']]
        closed = all(b['status'] == 'CONDITIONALLY_EXCLUDED_BRANCH' for b in branches)
        counts['CLOSED_CLASSES' if closed else 'OPEN_CLASSES'] += 1
        records.append(dict(cls=cls, status='CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES' if closed else 'OPEN', branches=branches))
    complete = selected == set(lookup) and all(r['status'] == 'CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES' for r in records) and not errors and not orphans
    return dict(status='CONDITIONALLY_COMPLETE_RESIDUAL_RANKONE' if complete else ('ERROR' if errors or orphans else 'OPEN'),
                classes=len(lookup), selected_classes=len(selected), counts=dict(counts), records=records,
                errors=errors, orphans=orphans, record_hashes=record_hashes, accepted_hashes=accepted,
                aggregation_source_sha256=sha(Path(__file__)),
                scope='Residual q1/q2 branches only; inherited classifications and other search branches remain separate dependencies.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--classes', type=int, nargs='+')
    ap.add_argument('--output', type=Path, default=HERE / 'aggregate.json')
    args = ap.parse_args()
    from rankone_reuse.pool import validate as validate_ready
    plan = validate_ready()
    out = aggregate(plan, args.classes)
    atomic_json(args.output, out)
    print(json.dumps({k: v for k, v in out.items() if k not in ('records', 'record_hashes', 'accepted_hashes')}), flush=True)


if __name__ == '__main__':
    main()
