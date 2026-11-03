"""Admit fresh GCC/Linux binaries with complete streams and genuine positives."""
import hashlib
import json
import platform
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from oracle_residual.common import HERE, source_hashes, activate, assignment, validate
from rankone_residual import common as base
from rankone_residual import controls
from rankone_residual.search import expected_slices
from mixed_frames.stage3.common import atomic_json, sha


def main():
    assert __debug__ and platform.system() == 'Linux'
    plan = validate()
    assignment()
    sources = {'compact': ROOT/'kernel_z_compact/z_enum_compact.c',
               'joint': ROOT/'kernel_z_joint/z_enum_joint.c'}
    commands = []
    for name, source in sources.items():
        cmd = ['cc', '-O3', '-std=c11', '-mcpu=native', str(source), '-o', str(HERE/f'linux_{name}')]
        subprocess.run(cmd, check=True)
        commands.append(cmd)
    activate()
    rows = {r['cls']: r for r in plan['records']}
    fixtures = json.loads((HERE/'fixtures.json').read_text())
    assert fixtures['status'] == 'COMPLETE_MAC_REFERENCE'
    for p, h in fixtures['reference_binaries'].items():
        assert sha(ROOT/p) == h
    records = []
    for f in fixtures['cases']:
        row = rows[f['cls']]
        s = next(x for x in expected_slices(plan, row, profiles=[f['profile']]) if x['oi'] == 0)
        model, data, H, order, cs, io = base.slice_setup(tuple(row['C']),
            tuple(s['branch']['N_orbits'][0]['N']), tuple(s['prefix']['H']),
            s['prefix']['p'], s['q'], s['t'])
        c = cs[s['ci']]
        ss = base.suffix_slice(model, s['branch']['N_orbits'][0]['N'], data, c)
        assert (H, order, s['ci'], base.digest(ss)) == (f['H'], f['order'], f['context'], f['orbit_digest'])
        assert base.payload(model, data, ss['first_representatives'][:64], c) == f['stdin']
        group = f'{len(c["group"])} 7\n' + '\n'.join(' '.join(map(str, a)) for a in c['group']) + '\n'
        assert group == f['group']
        with tempfile.TemporaryDirectory() as td:
            gp = Path(td)/'group.txt'
            gp.write_text(group)
            start = time.monotonic()
            p = subprocess.run([str(HERE/f'linux_{f["backend"]}'), *f['argv'],
                '--prefix-group', str(gp)], input=f['stdin'], text=True,
                capture_output=True, check=True, timeout=180)
            wall = time.monotonic()-start
        st = json.loads(p.stderr)
        assert {k:v for k,v in st.items() if k!='seconds'} == f['expected_stats']
        assert hashlib.sha256(p.stdout.encode()).hexdigest() == f['expected_stdout_sha256']
        r = dict(cls=f['cls'], profile=f['profile'], backend=f['backend'],
                 mac_wall=f['mac_wall'], linux_wall=wall, status='PASS', leaves=st['leaves'])
        records.append(r)
        print(json.dumps(r), flush=True)
    # The original control function uses both native backends and full 2^20
    # physical dual-word checks. Redirect its imported native entry point only.
    controls.z_search = base.z_search
    rng = random.Random(20261014)
    positives = []
    for hull in (1, 3):
        row = next(r for r in plan['records'] if r['hull'] == hull)
        for branch in row['branches']:
            for prefix in base.prefixes(plan, branch['q'], branch['t']):
                r = controls.control(row, branch, prefix, rng)
                assert r['status'] == 'PASS'
                positives.append(r)
                print(json.dumps(dict(stage='physical_positive', **r)), flush=True)
    assert len(positives) == 12
    report = dict(status='PASS_LINUX_NATIVE_AND_PHYSICAL_CONTROLS',
        source_hashes=source_hashes(), commands=commands, compiler=subprocess.check_output(['cc','--version'],text=True),
        runtime=dict(system=platform.platform(), python=sys.version), records=records,
        positives=positives, binary_hashes={str((HERE/f'linux_{n}').relative_to(ROOT)):
            sha(HERE/f'linux_{n}') for n in sources},
        scope='Portable recompilation admission, eight complete Mac stream/counter matches and twelve genuine distance-three positives. No distance-six target result.')
    atomic_json(HERE/'verification.json', report)
    print(json.dumps(dict(status=report['status'], comparisons=len(records), positives=len(positives))), flush=True)


if __name__ == '__main__':
    main()
