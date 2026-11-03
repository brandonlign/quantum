"""Independent complete small controls and physical positive fixed-Z controls for both new lanes.

Fixtures are the Stage-4 physical models (n=6 hull 0, n=7 hull 2 singular). For each mu:
  U = all leaves of the UNREDUCED frozen binary (canonical=0, domains=0) -- no symmetry at all;
  F = leaves of fast binary (canonical, frozen suffix-orbit root mask regenerated here, --group full closure).
Assert {joint_min(x) : x in U} == F exactly (complete orbit coverage, no duplicates), and also
check Stage-4-lane leaves S satisfy F subset S and joint_min(S)==F.
Positive controls: COMPLETE fixed-Z Y completion at every F leaf reproduces accepted physical solutions
whenever the unreduced lane has any. Union-budget lane: union-filter equality on fixtures with budgets.
"""
import itertools
import json
import time
from model import ExtensionModel, coordinates
from mixed_frames.frames import adapted_frame
from mixed_frames.stage3.common import input_text
from mixed_frames.stage4.symmetry import automorphisms, coverage
from mixed_frames.stage4.common import eligible
from mixed_frames.stage4.run import complete_leaf
from .common import FAST, FROZEN, HERE, TMP, group_closure, joint_min, run


def setup(m, f, mu, name):
    aut = automorphisms(m)
    maps = [g['syndrome_columns'] for g in aut['generators']]
    reps = {o['representative'] for o in coverage(eligible(m, f, mu), maps)}
    TMP.mkdir(exist_ok=True)
    mask = TMP / f'{name}.sym'
    mask.write_text(str(1 << m.b) + '\n' + ' '.join(str(int(a in reps)) for a in range(1 << m.b)) + '\n')
    elems = group_closure(maps, m.b)
    gpath = TMP / f'{name}.group'
    nontriv = [g for g in elems if list(g) != list(range(1 << m.b))]
    if nontriv:
        gpath.write_text(f'{1 << m.b} {len(nontriv)}\n' + '\n'.join(' '.join(map(str, g)) for g in nontriv) + '\n')
    else:
        gpath = None
    return mask, elems, gpath, aut


def inp(m, f, name, **kw):
    path = TMP / f'{name}.txt'
    path.write_text(input_text(m, f, **kw))
    return path


def main():
    report = dict(status='PASS', fixtures=[], positives=[], union=[])
    t = time.monotonic()
    for m in [ExtensionModel(6, 0, 3, 4, [15, 240], []), ExtensionModel(7, 1, 2, 4, [240, 15], [24, 3])]:
        f = adapted_frame(4, m.R)
        for mu in [0, m.d, m.d + 1]:
            name = f'fx_{m.n}_{mu}'
            mask, elems, gpath, aut = setup(m, f, mu, name)
            for rank in [False, True]:
                path = inp(m, f, name, mu=mu, rank=rank, domains=False)
                U = run(FROZEN, path, '--canonical', 0)
                S = run(FROZEN, path, '--canonical', 1, '--symmetry', mask)
                lanes = {}
                for partial in [0, 1]:
                    args = ['--canonical', 1, '--symmetry', mask]
                    if gpath:
                        args += ['--group', gpath, '--group-partial', partial]
                    lanes[partial] = run(FAST, path, *args)
                assert U['complete'] and S['complete'] and all(r['complete'] for r in lanes.values())
                mins = {tuple(joint_min(x, elems)) for x in U['leaves']}
                Sset = {tuple(x) for x in S['leaves']}
                for partial, F in lanes.items():
                    Fset = {tuple(x) for x in F['leaves']}
                    assert len(Fset) == len(F['leaves'])
                    assert Fset == mins, (name, rank, partial, sorted(mins ^ Fset))
                    assert Fset <= Sset
                assert {tuple(joint_min(x, elems)) for x in S['leaves']} == mins
                report['fixtures'].append(dict(n=m.n, hull=m.hull, mu=mu, rank=rank, group_order=len(elems),
                    unreduced_leaves=U['leaf_count'], unreduced_nodes=U['nodes'], stage4_leaves=S['leaf_count'],
                    stage4_nodes=S['nodes'], fast_leaves=lanes[1]['leaf_count'], fast_nodes=lanes[1]['nodes'],
                    fast_nodes_no_partial=lanes[0]['nodes'], orbits=len(mins)))
        # Physical positive: complete Y at every fast leaf (mu=d) vs. unreduced.
        mu = m.d; name = f'fx_{m.n}_{mu}'
        mask, elems, gpath, aut = setup(m, f, mu, name)
        path = inp(m, f, name, mu=mu, rank=True, domains=False)
        args = ['--canonical', 1, '--symmetry', mask] + (['--group', gpath] if gpath else [])
        F = run(FAST, path, *args); U = run(FROZEN, path, '--canonical', 0)
        u_ok = {}
        for x in U['leaves']:
            u_ok.setdefault(tuple(joint_min(x, elems)), []).append(len(complete_leaf(m, f, x, mu, None, m.d, m.d + 1)['accepted']))
        f_ok = {tuple(x): len(complete_leaf(m, f, x, mu, None, m.d, m.d + 1)['accepted']) for x in F['leaves']}
        assert set(f_ok) == set(u_ok)
        # A representative has a solution iff any orbit member does (physical invariance).
        for k, n in f_ok.items():
            assert (n > 0) == any(v > 0 for v in u_ok[k]), (k, n, u_ok[k])
        assert sum(f_ok.values()) > 0
        report['positives'].append(dict(n=m.n, hull=m.hull, fast_leaves=len(f_ok), fast_accepted=sum(f_ok.values()),
            unreduced_leaves=U['leaf_count'], unreduced_accepted=sum(sum(v) for v in u_ok.values()),
            orbits_with_solutions=sum(n > 0 for n in f_ok.values())))
        # Union-budget lane on fixture with exact profile budgets taken from actual measurements.
        meas = set()
        for x in F['leaves']:
            for s in complete_leaf(m, f, x, mu, None, m.d, m.d + 1)['measurements']:
                if s['minimum'] >= m.d:
                    meas.add((s['n6'], s['n7']))
        budgets = sorted(meas)[:4]
        if len(budgets) >= 2:
            L6 = max(b[0] for b in budgets); L7 = max(b[1] for b in budgets); L = max(sum(b) for b in budgets)
            upath = inp(m, f, name + '_u', mu=mu, rank=True, domains=False, counts=(L6, L7), w6=m.d, w7=m.d + 1)
            lines = upath.read_text().splitlines(); h = lines[5].split(); h[0] = str(L); lines[5] = ' '.join(h)
            upath.write_text('\n'.join(lines) + '\n')
            Uu = run(FAST, upath, '--relaxed', 1, '--canonical', 1, '--symmetry', mask)
            for bud in budgets:
                bpath = inp(m, f, name + f'_{bud[0]}{bud[1]}', mu=mu, rank=True, domains=False, counts=bud, w6=m.d, w7=m.d + 1)
                direct = run(FROZEN, bpath, '--canonical', 1, '--symmetry', mask)['leaves']
                keep = [x for x in Uu['leaves'] if run(FROZEN, bpath, '--prefix', ','.join(map(str, x)), '--canonical', 1, '--symmetry', mask)['leaves']]
                assert keep == direct, (bud, keep, direct)
                report['union'].append(dict(n=m.n, budget=bud, union=(L, L6, L7), union_leaves=Uu['leaf_count'],
                    direct_leaves=len(direct), union_nodes=Uu['nodes']))
    report['seconds'] = time.monotonic() - t
    (HERE / 'controls.json').write_text(json.dumps(report, indent=1) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
