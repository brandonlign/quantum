"""fast888_production: frozen preparation for the G+U superset census of 888 mu6 profiles.

Reads mixed_frames/stage4 and fast888_opus strictly read-only; writes only under fast888_production/.
"""
import json
from pathlib import Path
from mixed_frames.stage3.common import atomic_json, sha, target
from mixed_frames.stage4.common import hashes as stage4_hashes
from fast888_opus.common import group_closure

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
S4 = ROOT / 'mixed_frames/stage4'
OPUS = ROOT / 'fast888_opus'
BINARY = OPUS / 'enumerate_fast'
BINARY_SHA = 'f2f726455603f3ea1c4f447210a1d1a947dd5b076c740dc3680ba7c7f1264e2e'
PROFILES = ('10', '11', '20')            # all A6>=1 profiles; (0,1),(0,2) are mu7 legacy rows
LINE6 = {'10': '1 1 0 8 6 7 4', '11': '2 1 1 8 6 7 4', '20': '2 2 0 8 6 7 4'}
UNION_LINE6 = '2 2 1 8 6 7 4'            # componentwise max; fast binary needs --relaxed 1
ARGS = ['--relaxed', '1', '--canonical', '1', '--ack', '1']


def own_hashes():
    """Frozen fingerprints: old sources (stage4 manifest's own map), old manifest, new tool, own code."""
    old = json.loads((S4 / 'manifest.json').read_text())['hashes']
    tool = {str(p.relative_to(ROOT)): sha(p) for p in
            [OPUS / 'enumerate_fast', OPUS / 'enumerate_fast.cpp', OPUS / 'common.py', OPUS / 'controls.py',
             OPUS / 'controls.json', ROOT / 'mixed_frames/fixed_z.py', ROOT / 'mixed_frames/frames.py', ROOT / 'model.py']}
    own = {str(p.relative_to(ROOT)): sha(p) for p in sorted(HERE.glob('*.py'))}
    return dict(stage4_sources=old, stage4_manifest_sha256=sha(S4 / 'manifest.json'),
                stage4_exclusions_sha256=sha(S4 / 'exclusions.json'), tool=tool, own=own)


def union_text(index):
    """Assert the three frozen stage-4 inputs differ only in line 6, then substitute the union budget."""
    texts = {pf: (S4 / f'inputs/888_{index}_{pf}.txt').read_text().splitlines() for pf in PROFILES}
    base = texts['10']
    for pf, lines in texts.items():
        assert len(lines) == 7 and lines[:5] == base[:5] and lines[6] == base[6] and lines[5] == LINE6[pf], (index, pf)
        assert (S4 / f'inputs/888_{index}_{pf}.sym').read_bytes() == (S4 / f'inputs/888_{index}_10.sym').read_bytes()
    assert base[0].split()[3] == '6'
    return '\n'.join([*base[:5], UNION_LINE6, base[6]]) + '\n'


def group_text(index):
    aut = json.loads((S4 / f'certificates/aut_{index}.json').read_text())
    b = len(aut['P'])
    elems = group_closure([g['syndrome_columns'] for g in aut['generators']], b)
    assert len(elems) == aut['graph_group_order'], (index, len(elems))
    rest = [g for g in elems if list(g) != list(range(1 << b))]
    return (f'{1 << b} {len(rest)}\n' + '\n'.join(' '.join(map(str, g)) for g in rest) + '\n') if rest else None, len(elems)


def legacy_rows():
    rows = json.loads((S4 / 'manifest.json').read_text())['classes']
    return {r['key']: r for r in rows}


def checkpoint_status(key, first):
    p = S4 / 'runs' / key / str(first) / 'checkpoint.json'
    if not p.exists():
        return 'PENDING'
    try:
        return json.loads(p.read_text())['status']
    except (ValueError, KeyError):
        return 'UNREADABLE'
