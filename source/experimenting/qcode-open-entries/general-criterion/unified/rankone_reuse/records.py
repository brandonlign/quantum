"""Validate frozen legacy records and separately fingerprinted reuse records."""
import json
from pathlib import Path
from rankone_residual import common as baseline
from rankone_residual.aggregate import verify_chunk
from mixed_frames.stage3.common import sha

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LEGACY = ROOT / 'rankone_residual/runs'
RUNS = HERE / 'runs'


def pool_hashes():
    return {str(p.relative_to(ROOT)): sha(p) for p in
            [HERE / name for name in ('pool.py', 'records.py', 'aggregate.py',
                                      'verify_integration.py', 'integration_verification.json')]}


def accepted_hashes():
    from rankone_reuse.common import fingerprints
    return [baseline.fingerprints(), fingerprints()]


def verify_record(path, row, s, k, *, accepted=None, require_pool=True):
    from rankone_reuse.common import validate_reuse
    accepted = accepted_hashes() if accepted is None else accepted
    r = json.loads(Path(path).read_text())
    assert r['hashes'] in accepted, 'unknown producer fingerprints'
    legacy = r['hashes'] == accepted[0]
    if legacy:
        assert 'z_reuse' not in r and 'reuse_pool_hashes' not in r, 'legacy record with reuse metadata'
    elif require_pool:
        assert r.get('reuse_pool_hashes') == pool_hashes(), 'reuse pool fingerprints'
    status = verify_chunk(Path(path), row, s, k, r['hashes'])
    if 'z_reuse' in r:
        assert not legacy
        validate_reuse(Path(path), row, s, k, accepted_hashes=accepted)
    return status


def record_path(name, *, legacy=LEGACY, folder=RUNS):
    """Never mask a saved FOUND or a conflicting completed record."""
    old, new = Path(legacy) / name, Path(folder) / name
    if old.exists() and new.exists():
        r = json.loads(old.read_text())
        assert r['status'] not in ('FOUND', 'EXCLUDED_SLICE_CHUNK'), ('conflicting legacy record', name)
    return new if new.exists() else (old if old.exists() else None)
