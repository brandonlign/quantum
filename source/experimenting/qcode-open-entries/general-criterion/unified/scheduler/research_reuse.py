"""Versioned controller: retain eight-core allocation, replace residual producer."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scheduler import research_remaining as baseline
from mixed_frames.stage3.common import sha


def validate():
    baseline_validate()
    ready = json.loads((ROOT / 'rankone_reuse/READY.json').read_text())
    assert ready['status'] == 'READY_FOR_RESIDUAL_REUSE_SWEEP'
    for path, digest in ready['hashes'].items():
        assert sha(ROOT / path) == digest, path
    assert not (ROOT / 'scheduler/STOP_REQUESTED').exists()


baseline_validate = baseline.validate


def main():
    baseline.JOBS['residual'].update(module='rankone_reuse.pool', aggregate='rankone_reuse.aggregate',
        report='rankone_reuse/aggregate.json', log='rankone_reuse/pool.log', pid='rankone_reuse/pool.pid')
    baseline.validate = validate
    baseline.main()


if __name__ == '__main__':
    main()
