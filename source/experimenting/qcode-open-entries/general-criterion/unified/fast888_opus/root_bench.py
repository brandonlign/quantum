"""Per-root timing of the full-group lane on saved COMPLETE roots (bounded: --seconds cap per call).
Pending roots: optionally run BOTH frozen and fast lanes under the same cap (no live files touched).
Checks new leaves subset of saved/frozen leaves whenever both COMPLETE."""
import json
import random
import sys
from .common import FAST, FROZEN, HERE, S4, checkpoint, group_file, manifest, run

CAP = 15


def main():
    mode = sys.argv[1]; keys = sys.argv[2].split(','); n = int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    rng = random.Random(seed); rows = {r['key']: r for r in manifest()['classes']}; out = []
    for key in keys:
        row = rows[key]; i = row['class_index']; pf = ''.join(map(str, row['profile']))
        gpath, elems = group_file(i)
        roots = list(row['first_images']); rng.shuffle(roots)
        for f in roots[:n]:
            inp = S4 / row['input']; mask = S4 / row['symmetry']
            base = ['--prefix', f, '--canonical', 1, '--symmetry', mask]
            new = run(FAST, inp, *base, *(['--group', gpath] if gpath else []), '--seconds', CAP, timeout=CAP + 5)
            rec = dict(key=key, first=f, group_order=len(elems), new=dict(status=new['status'], nodes=new['nodes'], seconds=new['seconds'], leaves=new['leaf_count']))
            if mode == 'saved':
                v = checkpoint(i, pf, f)
                if v is None:
                    continue
                old = dict(status='COMPLETE', nodes=v['summary']['nodes'], seconds=v['summary']['seconds'], leaves=len(v['leaves']), X=v['leaves'])
            else:
                o = run(FROZEN, inp, *base, '--seconds', CAP, timeout=CAP + 5)
                old = dict(status=o['status'], nodes=o['nodes'], seconds=o['seconds'], leaves=o['leaf_count'], X=o['leaves'])
            if old['status'] == 'COMPLETE' and new['status'] == 'COMPLETE':
                assert set(map(tuple, new['leaves'])) <= set(map(tuple, old['X'])), (key, f)
            old.pop('X'); rec['old'] = old; out.append(rec); print(json.dumps(rec), flush=True)
    (HERE / f'root_bench_{mode}_{seed}_{"_".join(keys)}.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
