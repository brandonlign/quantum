"""Union-budget traversal for profiles (1,0),(1,1),(2,0): one Z tree, exact per-profile filter.

For each sampled root whose three live Stage-4 checkpoints are COMPLETE (read-only):
  * run the FROZEN Stage-4 binary once on the union input (L,L6,L7)=(3,2,1) [and optionally
    the fast binary with relaxed (2,2,1)];
  * replay every union leaf at full prefix under each profile's frozen input + mask;
  * assert the filtered set equals the saved checkpoint leaf list exactly, per profile.
Each subprocess has a hard timeout; roots above MAX_NODES saved are skipped (<=20 CPU s/call).
"""
import json
import random
import sys
import time
from .common import FAST, FROZEN, HERE, S4, PROFILE_LINE, TMP, checkpoint, manifest, run, union_input

MAX_NODES = 150_000


def relaxed_input(index):
    base = union_input(index).read_text().splitlines()
    base[5] = '2 2 1 8 6 7 4'
    path = TMP / f'relaxed_{index}.txt'
    path.write_text('\n'.join(base) + '\n')
    return path


def union_frozen_input(index):
    base = union_input(index).read_text().splitlines()
    base[5] = '3 2 1 8 6 7 4'
    path = TMP / f'unionfrozen_{index}.txt'
    path.write_text('\n'.join(base) + '\n')
    return path


def one(index, first, lanes=('frozen', 'relaxed')):
    saved = {pf: checkpoint(index, pf, first) for pf in PROFILE_LINE}
    if any(v is None for v in saved.values()):
        return None
    mask = S4 / f'inputs/888_{index}_10.sym'
    out = dict(class_index=index, first=first,
               saved={pf: dict(nodes=v['summary']['nodes'], leaves=len(v['leaves']), seconds=v['summary']['seconds'])
                      for pf, v in saved.items()})
    for lane in lanes:
        if lane == 'frozen':
            binary, inp = FROZEN, union_frozen_input(index)
            r = run(binary, inp, '--prefix', first, '--canonical', 1, '--symmetry', mask)
        else:
            binary, inp = FAST, relaxed_input(index)
            r = run(binary, inp, '--relaxed', 1, '--prefix', first, '--canonical', 1, '--symmetry', mask)
        assert r['complete'], r
        union = [tuple(x) for x in r['leaves']]
        assert len(set(union)) == len(union)
        filt_seconds = 0.0
        for pf, v in saved.items():
            keep = []
            for x in union:
                t = run(FROZEN, S4 / f'inputs/888_{index}_{pf}.txt', '--prefix', ','.join(map(str, x)),
                        '--canonical', 1, '--symmetry', mask)
                assert t['complete']
                filt_seconds += t['seconds']
                assert t['leaves'] in ([], [list(x)]), t
                if t['leaves']:
                    keep.append(list(x))
            # Exact equality with the live COMPLETE checkpoint, order-preserving (DFS order).
            assert keep == v['leaves'], (index, first, pf, lane, keep, v['leaves'])
        out[lane] = dict(nodes=r['nodes'], leaves=len(union), seconds=r['seconds'],
                         candidates=r['candidates'], filter_binary_seconds=filt_seconds)
    return out


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    per_class = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    rng = random.Random(seed)
    cands = json.loads((HERE / 'candidates.json').read_text())
    picks = [(c[0], c[1]) for c in cands['leafy'] if max(v[0] for v in c[2].values()) <= MAX_NODES]
    man = manifest()
    for row in man['classes']:
        if row['profile'] != [1, 0]:
            continue
        roots = list(row['first_images'])
        rng.shuffle(roots)
        got = 0
        for f in roots:
            p = S4 / 'runs' / row['key'] / str(f) / 'checkpoint.json'
            if not p.exists():
                continue
            c = json.loads(p.read_text())
            if c['status'] == 'COMPLETE' and c['enumerator']['nodes'] <= MAX_NODES:
                picks.append((row['class_index'], f)); got += 1
            if got >= per_class:
                break
    picks = sorted(set(picks))
    started = time.monotonic(); rows = []
    for i, f in picks:
        r = one(i, f)
        if r is None:
            continue
        rows.append(r)
        print(json.dumps(r), flush=True)
    tot = dict(roots=len(rows),
               saved_nodes=sum(sum(v['nodes'] for v in r['saved'].values()) for r in rows),
               saved_seconds=sum(sum(v['seconds'] for v in r['saved'].values()) for r in rows),
               saved_leaves=sum(sum(v['leaves'] for v in r['saved'].values()) for r in rows))
    for lane in ('frozen', 'relaxed'):
        tot[lane] = {k: sum(r[lane][k] for r in rows) for k in ('nodes', 'leaves', 'seconds', 'filter_binary_seconds')}
    tot['wall'] = time.monotonic() - started
    (HERE / f'union_bench_{seed}.json').write_text(json.dumps(dict(totals=tot, rows=rows), indent=1) + '\n')
    print(json.dumps(tot), flush=True)


if __name__ == '__main__':
    main()
