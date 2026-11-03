"""Strict read-only re-aggregation of mixed_frames/stage2 evidence (root session).

Closes the three gaps reported in sol_workers/E_stage2_audit.md without touching
the running pipeline: (1) the enumerator log must end in exactly one COMPLETE
summary equal to the checkpoint's; (2) surviving solutions are recomputed from
the completion files, not trusted from the checkpoint; (3) every completion's
mu_out must equal the class's mu_out. A class counts as excluded only if every
manifest first image has a partition passing all checks with zero survivors.
"""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]/'stage2'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_partition(folder, row):
    cp = json.loads((folder/'checkpoint.json').read_text())
    if cp.get('status') != 'COMPLETE':
        return False, 'not complete', 0
    lines = [json.loads(s) for s in (folder/'enumerator.jsonl').read_text().splitlines() if s.strip()]
    summaries = [o for o in lines if o['type'] == 'summary']
    leaves = [o['X'] for o in lines if o['type'] == 'leaf']
    if len(summaries) != 1 or lines[-1]['type'] != 'summary':
        return False, 'summary missing or not last', 0
    s = summaries[0]
    if not (s['status'] == 'COMPLETE' and s['complete'] and s['leaf_count'] == len(leaves)):
        return False, 'summary not COMPLETE', 0
    if s != cp['enumerator'] or sha(folder/'enumerator.jsonl') != cp['log_sha256']:
        return False, 'summary/log mismatch with checkpoint', 0
    if len({tuple(x) for x in leaves}) != len(leaves):
        return False, 'duplicate leaves', 0
    survivors = 0
    for xs in leaves:
        key = hashlib.sha256(json.dumps(xs).encode()).hexdigest()
        f = folder/'completions'/(key+'.json')
        if not f.exists():
            return False, 'missing completion', 0
        c = json.loads(f.read_text())
        if c['X'] != xs or c['mu_out'] != row['mu_out'] or not c['complete'] or c['status'] != 'COMPLETE':
            return False, 'bad completion record', 0
        if not (c['completion']['complete'] and c['completion']['status'] == 'COMPLETE'):
            return False, 'inner completion incomplete', 0
        survivors += len(c['accepted'])
    return True, 'ok', survivors


def main():
    manifest = json.loads((HERE/'manifest.json').read_text())
    out, summary = [], {}
    for row in manifest['classes']:
        good, surv, bad = 0, 0, []
        for first in row['first_images']:
            folder = HERE/'runs'/f"{row['frame']}_{row['class_index']}"/str(first)
            if not (folder/'checkpoint.json').exists():
                bad.append((first, 'pending')); continue
            ok, why, s = check_partition(folder, row)
            if ok:
                good += 1; surv += s
            else:
                bad.append((first, why))
        complete = good == len(row['first_images'])
        excluded = complete and surv == 0
        out.append(dict(frame=row['frame'], class_index=row['class_index'], mu_out=row['mu_out'],
                        partitions=len(row['first_images']), verified_partitions=good,
                        survivors=surv, complete=complete, excluded=excluded,
                        problems=[b for b in bad if b[1] != 'pending' and b[1] != 'not complete'][:5]))
        t = summary.setdefault(row['frame'], dict(classes=0, excluded=0, found=0))
        t['classes'] += 1; t['excluded'] += excluded; t['found'] += bool(surv)
    res = dict(status='STRICT_REAGGREGATION', manifest_sha256=sha(HERE/'manifest.json'),
               summary=summary, classes=out,
               scope='Mixed-frame (668/778 mu_out=8, 688 mu=6) recentred searches only; '
                     'conditional on the 56-class length-nine census; not a [[17,3,6]] result.')
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('strict_aggregate.json')
    target.write_text(json.dumps(res, indent=1))
    print(json.dumps(summary))
    print('problems:', [(c['frame'], c['class_index'], c['problems']) for c in out if c['problems']])


if __name__ == '__main__':
    main()
