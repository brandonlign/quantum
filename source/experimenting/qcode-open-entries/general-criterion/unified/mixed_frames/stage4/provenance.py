"""Replay frozen legacy exclusion evidence, redirecting all outputs to Stage 4."""
import json
from .common import HERE,atomic_json,sha
from mixed_frames.stage3 import overlap


def main():
    source=HERE.parent/'stage3/suffix_tables.json'
    (HERE/'suffix_tables.json').write_bytes(source.read_bytes())
    original=overlap.HERE
    try:
        overlap.HERE=HERE
        overlap.main()
    finally:overlap.HERE=original
    assert json.loads((HERE/'overlap.json').read_text())==json.loads((HERE.parent/'stage3/overlap.json').read_text())
    atomic_json(HERE/'provenance.json',dict(status='PASS',
        frozen_stage3_overlap_sha256=sha(HERE.parent/'stage3/overlap.json'),
        replay_sha256=sha(HERE/'overlap.json'),suffix_tables_sha256=sha(source),
        scope='Independent Stage-3 evidence replay with output relocated to Stage 4; no legacy DFS rerun.'))


if __name__=='__main__':main()
