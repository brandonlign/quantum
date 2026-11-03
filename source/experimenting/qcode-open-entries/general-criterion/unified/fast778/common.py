"""Isolated frozen Stage-2 778 inputs using the validated S6 backend."""
import json
from pathlib import Path
from mixed_frames.stage3.common import atomic_json,budget,graph_counts,sha,target
from mixed_frames.stage2.common import hashes as original_hashes,input_text

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
BINARY=ROOT/'speedup_review/enumerate_688'


def case_input(model,frame,profile):
    assert profile==[0,2] or profile==(0,2)
    return input_text(model,frame,8)


def hashes():
    out=original_hashes()
    paths=[HERE/'common.py',HERE/'run.py',HERE/'verify_pipeline.py',HERE/'verify.py',
           HERE/'PROOF.md',BINARY,ROOT/'speedup_review/enumerate_688.cpp',
           ROOT/'mixed_frames/stage3/common.py',ROOT/'mixed_frames/stage3/run.py',
           ROOT/'mixed_frames/stage3/frames.py',ROOT/'speedup_review/PROOF.md',
           ROOT/'speedup_review/verify_688_symmetry.py']
    for p in paths:out[str(p.relative_to(ROOT))]=sha(p)
    return out


def prepare():
    saved=json.loads((ROOT/'mixed_frames/stage2/manifest.json').read_text())
    assert saved['hashes']==original_hashes()
    (HERE/'inputs').mkdir(exist_ok=True)
    rows=[]
    for old in saved['classes']:
        if old['frame']!='778':continue
        row=dict(old,key=f'778_{old["class_index"]}_02',profile=[0,2],budget=[0,0])
        data=(ROOT/'mixed_frames/stage2'/row['input']).read_text()
        assert sha(ROOT/'mixed_frames/stage2'/row['input'])==row['input_sha256']
        (HERE/row['input']).write_text(data)
        rows.append(row)
    assert len(rows)==52 and sum(len(r['first_images']) for r in rows)==650
    atomic_json(HERE/'manifest.json',dict(classes=rows,hashes=hashes(),partitions=650,
        source_manifest_sha256=sha(ROOT/'mixed_frames/stage2/manifest.json')))
    return dict(cases=52,roots=650)


def validate_manifest():
    saved=json.loads((HERE/'manifest.json').read_text())
    assert saved['hashes']==hashes()
    assert saved['source_manifest_sha256']==sha(ROOT/'mixed_frames/stage2/manifest.json')
    source=json.loads((ROOT/'mixed_frames/stage2/manifest.json').read_text())
    required=[dict(old,key=f'778_{old["class_index"]}_02',profile=[0,2],budget=[0,0]) for old in source['classes'] if old['frame']=='778']
    assert saved['classes']==required and saved['partitions']==650
    for name in ('verification.json','pipeline_verification.json'):
        v=json.loads((HERE/name).read_text());assert v['status']=='PASS' and v['hashes']==saved['hashes']
    for r in saved['classes']:
        m,f=target('778',r['class_index'])
        assert (HERE/r['input']).read_text()==case_input(m,f,r['profile'])
        assert sha(HERE/r['input'])==r['input_sha256']
        assert r['budget']==list(budget(f,r['profile']))==[0,0]
    return saved


if __name__=='__main__':print(json.dumps(prepare()))
