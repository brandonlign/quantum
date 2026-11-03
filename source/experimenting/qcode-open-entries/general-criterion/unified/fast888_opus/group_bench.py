"""Full suffix-group x prefix-S8 joint-lex canonical lane versus live COMPLETE Stage-4 checkpoints.

Whole-case test: every root of a COMPLETE case is rerun with --group (same frozen input and mask).
Checks: (a) new leaves are a subset of saved leaves of that root (only extra pruning);
(b) leaf-orbit coverage: the set of joint minima of all saved leaves of the case equals the set
of new leaves (each new leaf is its own joint min, each saved orbit hit); (c) zero-leaf roots
match. Every subprocess <= 20 s; cases with a saved root above MAX_NODES are skipped.
"""
import json
import sys
import time
from .common import FAST, S4, HERE, checkpoint, group_file, joint_min, manifest, run

MAX_NODES = 1_200_000


def case(row, partial=1):
    i = row['class_index']; pf = ''.join(map(str, row['profile']))
    gpath, elems = group_file(i)
    saved = {f: checkpoint(i, pf, f) for f in row['first_images']}
    if any(v is None for v in saved.values()):
        return dict(key=row['key'], skipped='incomplete')
    if max(v['summary']['nodes'] for v in saved.values()) > MAX_NODES:
        return dict(key=row['key'], skipped='large root')
    inp = S4 / row['input']; mask = S4 / row['symmetry']
    old_leaves = []; new_leaves = []; old_nodes = new_nodes = 0; old_s = new_s = 0.0
    zero_match = 0
    for f, v in saved.items():
        args = ['--prefix', f, '--canonical', 1, '--symmetry', mask]
        if gpath is not None:
            args += ['--group', gpath, '--group-partial', partial]
        r = run(FAST, inp, *args, timeout=20)
        assert r['complete'], r
        assert set(map(tuple, r['leaves'])) <= set(map(tuple, v['leaves'])), (row['key'], f)
        zero_match += (not r['leaves']) == (not v['leaves'])
        old_leaves += v['leaves']; new_leaves += r['leaves']
        old_nodes += v['summary']['nodes']; new_nodes += r['nodes']
        old_s += v['summary']['seconds']; new_s += r['seconds']
    mins = {tuple(joint_min(x, elems)) for x in old_leaves}
    new = {tuple(x) for x in new_leaves}
    assert len(new) == len(new_leaves)
    assert new == mins, (row['key'], sorted(mins - new), sorted(new - mins))
    return dict(key=row['key'], group_order=len(elems), roots=len(saved), saved_nodes=old_nodes, new_nodes=new_nodes,
                saved_seconds=round(old_s, 3), new_seconds=round(new_s, 3), saved_leaves=len(old_leaves),
                new_leaves=len(new_leaves), leaf_orbits=len(mins), zero_leaf_roots_with_zero=zero_match)


def main():
    keys = sys.argv[1:]
    man = manifest(); out = []
    t = time.monotonic()
    for row in man['classes']:
        if row['key'] in keys:
            r = case(row); out.append(r); print(json.dumps(r), flush=True)
    (HERE / ('group_bench_' + '_'.join(keys) + '.json')).write_text(json.dumps(dict(rows=out, wall=time.monotonic() - t), indent=1) + '\n')


if __name__ == '__main__':
    main()
