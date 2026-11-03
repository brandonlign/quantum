"""Bounded exact native/reused streams, own-H Y evidence, and failure paths."""
import copy
import json
import shutil
import sys
import time
from pathlib import Path
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / 'm8r1')]
from rankone_residual import common as base
from rankone_residual.search import expected_slices
from rankone_reuse import common
from mixed_frames.stage3.common import atomic_json, sha


def scrub(x):
    if isinstance(x, dict):
        return {k: scrub(v) for k, v in x.items() if k != 'seconds'}
    if isinstance(x, list):
        return [scrub(v) for v in x]
    return x


def main():
    start = time.monotonic()
    plan = base.validate(require_controls=True)
    lookup = {r['cls']: r for r in plan['records']}
    hashes = common.source_hashes()
    folder = HERE / 'engine_control'
    assert not folder.exists(), 'Preserve prior verification; create a new version after changes'
    folder.mkdir()
    pairs = []
    for cls, k in [(0, 0), (0, 1), (1, 0), (1, 3)]:
        row = lookup[cls]
        s = next(expected_slices(plan, row, profiles=[2]))
        source = ROOT / 'rankone_residual/runs' / base.task_name(cls, 1, 0, 0, 1, 0, k)
        saved = json.loads(source.read_text())
        reused = common.run_task(row, s, k, hashes, source_path=source)
        native = common.run_task(row, s, k, hashes)
        for field in ('z', 'leaves', 'tally', 'H', 'first', 'status', 'backend'):
            assert scrub(reused[field]) == scrub(native[field]), (cls, k, field)
        assert reused['status'] == native['status'] == 'EXCLUDED_SLICE_CHUNK'
        path = folder / s['names'][k]
        atomic_json(path, reused)
        common.validate_reuse(path, row, s, k, [base.fingerprints(), hashes])
        pairs.append(dict(cls=cls, chunk=k, leaves=len(reused['leaves']),
                          native_seconds=native['seconds'], reused_seconds=reused['seconds'],
                          source_seconds=saved['seconds'], source_sha256=sha(source),
                          result_sha256=sha(path)))
    row = lookup[0]
    s = next(expected_slices(plan, row, profiles=[2]))
    source = ROOT / 'rankone_residual/runs' / base.task_name(0, 1, 0, 0, 1, 0, 0)
    original_z = base.z_search
    rejected = []
    for kind in ('different input', 'exception during Y'):
        try:
            if kind == 'different input':
                with patch.object(common, 'native_input', side_effect=[dict(a=1), dict(a=2)]):
                    common.run_task(row, s, 0, hashes, source_path=source)
            else:
                with patch.object(base, 'solve', side_effect=RuntimeError('injected Y failure')):
                    common.run_task(row, s, 0, hashes, source_path=source)
        except (AssertionError, RuntimeError):
            rejected.append(kind)
        else:
            raise AssertionError(('failure not rejected', kind))
        assert base.z_search is original_z
    # Test the unchanged raw-candidate path using an explicitly injected FOUND.
    # This is a control-flow test, not an actual distance-six construction.
    fake = dict(status='FOUND', reason=None, solutions=[dict(model_columns=[0] * 16)])
    candidates = folder / 'candidates'
    with patch.object(base, 'solve', return_value=copy.deepcopy(fake)), \
         patch.object(base, 'physical_witness', side_effect=RuntimeError('injected physical verifier failure')):
        found = common.run_task(row, s, 0, hashes, source_path=source, candidate_dir=candidates)
    assert found['status'] == 'FOUND' and 'physical_error' in found
    raw = json.loads((ROOT / found['candidate_file']).read_text())
    assert raw['record']['status'] == 'FOUND' and raw['record']['z_reuse'] == found['z_reuse']
    assert base.z_search is original_z
    report = dict(status='PASS', pairs=pairs, rejected=rejected,
                  own_H_Y_evidence_matches_native=True, restored_after_failure=True,
                  injected_FOUND_persisted_before_failed_physical_check=True,
                  positive_control_scope='injected FOUND tests persistence only; baseline genuine d=3 physical controls remain inherited',
                  hashes=common.source_hashes(), seconds=time.monotonic() - start)
    atomic_json(HERE / 'engine_verification.json', report)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
