"""Separately admitted Linux producer; original binaries remain reference bytes."""
import json
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'vendor'))
from rankone_residual import common as base
from mixed_frames.stage3.common import sha

ROOT = HERE.parent
FILES = ['oracle_residual/common.py', 'oracle_residual/prepare.py',
         'oracle_residual/verify.py', 'oracle_residual/pool.py',
         'oracle_residual/aggregate.py', 'oracle_residual/bridge.py',
         'oracle_residual/launch.py', 'rankone_residual/controls.py',
         'fixed_kernel/verify_large.py', 'gauge_fix.py',
         'rankone_residual/READY.json', 'rankone_residual/controls.json',
         'oracle_residual/vendor/full.py', 'oracle_residual/vendor/sympl.py',
         'oracle_residual/vendor/xstage.py', 'oracle_residual/vendor/zstage.py',
         'oracle_residual/vendor/release16_cells.json', 'oracle_residual/external_paths.json',
         'oracle_residual/PROOF.md', 'scheduler/research_oracle.py',
         'oracle_residual/assignment.json', 'oracle_residual/fixtures.json']


def source_hashes():
    return dict(base.fingerprints(), **{p: sha(ROOT / p) for p in FILES})


def fingerprints():
    report = HERE / 'verification.json'
    d = json.loads(report.read_text())
    assert d['status'] == 'PASS_LINUX_NATIVE_AND_PHYSICAL_CONTROLS'
    assert d['source_hashes'] == source_hashes(), 'Linux admission source changes'
    for p, h in d['binary_hashes'].items():
        assert sha(ROOT / p) == h, ('Linux binary', p)
    return dict(source_hashes(), **d['binary_hashes'],
                **{'oracle_residual/verification.json': sha(report)})


def activate():
    base.Z_BINARY = HERE / 'linux_compact'
    base.Z_JOINT = HERE / 'linux_joint'


def validate():
    """Relocate only explicitly listed external hash inputs, preserving bytes."""
    paths = json.loads((HERE / 'external_paths.json').read_text())
    saved = base.sha

    def portable_sha(path):
        entry = paths.get(str(path))
        if entry is None:
            return sha(path)
        target = ROOT / entry['copy']
        digest = sha(target)
        assert digest == entry['sha256'], 'external dependency bytes changed'
        return digest

    try:
        base.sha = portable_sha
        return base.validate(require_controls=True)
    finally:
        base.sha = saved


def assignment():
    d = json.loads((HERE / 'assignment.json').read_text())
    assert d['status'] == 'DISJOINT_CLASS_PARTITION'
    plan = json.loads(base.PLAN.read_text())
    universe = {r['cls'] for r in plan['records']}
    local, remote = set(d['local']), set(d['remote'])
    assert not local & remote and local | remote == universe
    assert len(local) == len(d['local']) and len(remote) == len(d['remote'])
    return d
