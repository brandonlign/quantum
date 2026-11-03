"""Frozen isolated inputs and implementation fingerprints for the S6 search."""
import json
from pathlib import Path
from mixed_frames.stage3.common import (
    atomic_json, budget, case_input, graph_counts, sha, target,
    hashes as original_hashes)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BINARY = HERE/'enumerate_688'


def hashes():
    out = original_hashes()
    for name in ('common.py','run.py','verify_backend.py','verify_688_symmetry.py',
                 'enumerate_688.cpp','enumerate_688','PROOF.md','verify_pipeline.py'):
        p = HERE/name
        out[str(p.relative_to(ROOT))] = sha(p)
    return out


def prepare():
    saved = json.loads((ROOT/'mixed_frames/stage3/manifest.json').read_text())
    assert original_hashes() == saved['hashes']
    verification = json.loads((HERE/'backend_verification.json').read_text())
    assert verification['status']=='PASS'
    assert verification['backend_sha256']==sha(BINARY)
    assert verification['source_sha256']==sha(HERE/'enumerate_688.cpp')
    (HERE/'inputs').mkdir(exist_ok=True)
    rows = []
    for old in saved['classes']:
        if old['frame']!='688':continue
        row = dict(old)
        data = (ROOT/'mixed_frames/stage3'/row['input']).read_text()
        assert sha(ROOT/'mixed_frames/stage3'/row['input'])==row['input_sha256']
        (HERE/row['input']).write_text(data)
        rows.append(row)
    assert len(rows)==156 and sum(len(r['first_images']) for r in rows)==1950
    atomic_json(HERE/'manifest.json',dict(classes=rows,hashes=hashes(),
        source_manifest_sha256=sha(ROOT/'mixed_frames/stage3/manifest.json'),
        partitions=1950,scope='S6-reduced 688; every certified first image retained'))
    return dict(cases=len(rows),roots=1950)


def validate_manifest():
    saved = json.loads((HERE/'manifest.json').read_text())
    assert saved['hashes']==hashes(), 'Frozen implementation changed'
    assert saved['source_manifest_sha256']==sha(ROOT/'mixed_frames/stage3/manifest.json')
    original = json.loads((ROOT/'mixed_frames/stage3/manifest.json').read_text())
    assert saved['classes']==[r for r in original['classes'] if r['frame']=='688']
    assert saved['partitions']==sum(len(r['first_images']) for r in saved['classes'])==1950
    backend=json.loads((HERE/'backend_verification.json').read_text())
    assert backend['status']=='PASS' and backend['backend_sha256']==sha(BINARY)
    assert backend['source_sha256']==sha(HERE/'enumerate_688.cpp')
    pipeline=json.loads((HERE/'pipeline_verification.json').read_text())
    assert pipeline['status']=='PASS' and pipeline['hashes']==saved['hashes']
    for r in saved['classes']:
        assert sha(HERE/r['input'])==r['input_sha256']
        m,f=target(r['frame'],r['class_index'])
        assert (HERE/r['input']).read_text()==case_input(m,f,r['profile'])
        assert r['budget']==list(budget(f,r['profile']))
    return saved


if __name__=='__main__':print(json.dumps(prepare()))
