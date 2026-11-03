"""Construct portable, exact-byte native fixtures on the running Mac."""
import hashlib
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from rankone_residual import common as base
from rankone_residual.search import expected_slices
from mixed_frames.stage3.common import atomic_json, sha


def main():
    plan = base.validate(require_controls=True)
    ids = sorted(r['cls'] for r in plan['records'])
    assignment = ROOT / 'oracle_residual/assignment.json'
    assert not assignment.exists(), 'Assignment is immutable once prepared'
    atomic_json(assignment, dict(status='DISJOINT_CLASS_PARTITION',
                local=ids[:-4096], remote=list(reversed(ids[-4096:]))))
    rows = {r['cls']: r for r in plan['records']}
    cases = []
    for cls, profile in [(0, p) for p in (0, 1, 3, 8, 12)] + [(5, p) for p in (4, 5, 13)]:
        row = rows[cls]
        s = next(x for x in expected_slices(plan, row, profiles=[profile]) if x['oi'] == 0)
        model, data, H, order, cs, io = base.slice_setup(tuple(row['C']),
            tuple(s['branch']['N_orbits'][0]['N']), tuple(s['prefix']['H']),
            s['prefix']['p'], s['q'], s['t'])
        c = cs[s['ci']]
        ss = base.suffix_slice(model, s['branch']['N_orbits'][0]['N'], data, c)
        reps = ss['first_representatives'][:64]
        stdin = base.payload(model, data, reps, c)
        group = f'{len(c["group"])} 7\n' + '\n'.join(' '.join(map(str, a)) for a in c['group']) + '\n'
        backend = 'joint' if profile in (0, 8) else 'compact'
        argv = ['--gamma', '--avoid-radical', '--first-orbits', '--kernel-data']
        if backend == 'joint':
            argv += ['--mitm4', '65536']
        with tempfile.TemporaryDirectory() as td:
            gp = Path(td) / 'group.txt'
            gp.write_text(group)
            started = time.monotonic()
            p = subprocess.run([str(base.Z_JOINT if backend == 'joint' else base.Z_BINARY),
                *argv, '--prefix-group', str(gp)], input=stdin, text=True,
                capture_output=True, check=True, timeout=120)
            wall = time.monotonic() - started
        st = json.loads(p.stderr)
        assert st['status'] == 'COMPLETE' and st.get('reason') is None
        cases.append(dict(cls=cls, profile=profile, backend=backend, argv=argv,
            stdin=stdin, group=group, expected_stdout_sha256=hashlib.sha256(p.stdout.encode()).hexdigest(),
            expected_stats={k: v for k, v in st.items() if k != 'seconds'},
            mac_wall=wall, orbit_digest=base.digest(ss), H=H, order=order,
            context=s['ci']))
    atomic_json(ROOT / 'oracle_residual/fixtures.json', dict(status='COMPLETE_MAC_REFERENCE',
        cases=cases, reference_binaries={str(p.relative_to(ROOT)): sha(p)
            for p in (base.Z_BINARY, base.Z_JOINT)},
        scope='Eight complete native streams and finite-integer setup identities; no new target exclusion'))
    print(json.dumps(dict(fixtures=len(cases), local=len(ids)-4096, remote=4096)), flush=True)


if __name__ == '__main__':
    main()
