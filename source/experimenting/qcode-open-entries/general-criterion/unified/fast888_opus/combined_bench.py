"""Combined lane (full suffix group + union budget) on a whole class triple 10/11/20.
For each profile: {replayed union leaves passing that profile} == {joint_min(saved leaves of that profile)}."""
import json
import sys
import time
from .common import FAST, HERE, S4, PROFILE_LINE, checkpoint, group_file, joint_min, manifest, run
from .union_bench import relaxed_input


def main():
    i = int(sys.argv[1]); t = time.monotonic()
    gpath, elems = group_file(i)
    row = next(r for r in manifest()['classes'] if r['key'] == f'888_{i}_10')
    mask = S4 / row['symmetry']; g = ['--group', gpath] if gpath else []
    saved_nodes = saved_s = new_nodes = new_s = 0; union = []; saved = {pf: [] for pf in PROFILE_LINE}
    for f in row['first_images']:
        for pf in PROFILE_LINE:
            v = checkpoint(i, pf, f); assert v is not None
            saved[pf] += v['leaves']; saved_nodes += v['summary']['nodes']; saved_s += v['summary']['seconds']
        r = run(FAST, relaxed_input(i), '--relaxed', 1, '--prefix', f, '--canonical', 1, '--symmetry', mask, *g, timeout=20)
        assert r['complete']; union += r['leaves']; new_nodes += r['nodes']; new_s += r['seconds']
    out = dict(class_index=i, group_order=len(elems), roots=len(row['first_images']), saved_nodes_3profiles=saved_nodes,
               saved_seconds_3profiles=round(saved_s, 3), combined_nodes=new_nodes, combined_seconds=round(new_s, 3),
               union_leaves=len(union), profiles={})
    for pf in PROFILE_LINE:
        keep = {tuple(x) for x in union if run(FAST, S4 / f'inputs/888_{i}_{pf}.txt', '--prefix', ','.join(map(str, x)),
                                              '--canonical', 1, '--symmetry', mask, *g)['leaves']}
        mins = {tuple(joint_min(x, elems)) for x in saved[pf]}
        assert keep == mins, (pf, sorted(keep ^ mins))
        out['profiles'][pf] = dict(saved_leaves=len(saved[pf]), orbits=len(mins))
    out['wall'] = time.monotonic() - t
    (HERE / f'combined_{i}.json').write_text(json.dumps(out, indent=1) + '\n'); print(json.dumps(out))


if __name__ == '__main__':
    main()
