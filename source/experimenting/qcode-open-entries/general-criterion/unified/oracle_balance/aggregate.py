"""Whole-class union, including optional Mac assistance on the Oracle shard."""
import json
import sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from oracle_residual.common import fingerprints, assignment
from oracle_residual.aggregate import verify
from rankone_reuse import aggregate as original
from rankone_reuse.records import accepted_hashes
from rankone_reuse.pool import validate
from mixed_frames.stage3.common import atomic_json, sha


def aggregate(selected=None):
    plan=validate()
    a=assignment()
    universe=set(a['local']+a['remote'])
    selected=universe if selected is None else set(selected)
    assert selected<=universe
    known=accepted_hashes()+[fingerprints()]
    # Each producer is strictly checked on its own. Only whole-class closure
    # from either producer is used; partial chunks are never silently combined.
    mac=original.aggregate(plan,selected)
    saved=original.verify_record
    original.verify_record=verify
    try:
        vm=original.aggregate(plan,set(a['remote']) & selected,legacy=ROOT/'oracle_residual/empty_legacy',
            folder=ROOT/'oracle_residual/runs',accepted=known)
    finally:
        original.verify_record=saved
    lookup={r['cls']:r for r in vm['records']}
    closed='CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES'
    records=[]
    for r in mac['records']:
        v=lookup.get(r['cls'])
        records.append(dict(r if r['status']==closed or v is None or v['status']!=closed else v,
            closure_producer='mac' if r['status']==closed else
                ('oracle' if v is not None and v['status']==closed else None)))
    errors=mac['errors']+vm['errors']
    for directory in ['rankone_residual/candidates','rankone_reuse/runs/candidates','oracle_residual/candidates']:
        for path in (ROOT/directory).glob('*.json'):
            errors.append(dict(path=str(path),error='Persisted raw candidate requires inspection'))
    orphans=mac['orphans']+vm['orphans']
    complete=selected==universe and len(records)==len(selected) and all(r['status']==closed for r in records) and not errors and not orphans
    out=dict(status='CONDITIONALLY_COMPLETE_RESIDUAL_RANKONE' if complete else
        ('ERROR' if errors or orphans else 'OPEN'), classes=len(plan['records']),selected_classes=len(selected),
        counts=dict(Counter(r['status'] for r in records)),records=records,
        errors=errors,orphans=orphans,
        record_hashes=dict(mac['record_hashes'],**vm['record_hashes']),
        accepted_hashes=known,source_hashes={str(p.relative_to(ROOT)):sha(p)
            for p in [Path(__file__),ROOT/'oracle_balance/watch.py',ROOT/'oracle_balance/stop_remote.py']},
        scope='Strict whole-class union of Mac and Oracle residual branches. Both producers are audited even for duplicate classes; a FOUND/error anywhere prevents closure. Other branches and final theorem review remain separate.')
    atomic_json(ROOT/'oracle_balance/aggregate.json',out)
    return out


if __name__=='__main__':
    r=aggregate()
    print(json.dumps({k:v for k,v in r.items() if k not in ('records','record_hashes','accepted_hashes')}),flush=True)
