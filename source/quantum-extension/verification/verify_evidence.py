#!/usr/bin/env python3
"""Bind the unified manuscript to frozen proof evidence without changing it.

Default: replay exact low-weight certificates; verify classification provenance,
full saved branch coverage and current implementation fingerprints; replay all
mixed-frame saved partitions; inspect residual coverage and a fresh sample.
--deep additionally repeats the full 2.4-million-record residual aggregation.
This checks evidence and arithmetic, not the truth of handwritten arguments.
"""
from __future__ import annotations
import argparse
import contextlib
import hashlib
import importlib
import io
import json
import pathlib
import subprocess
import sys
import tempfile
import time
from collections import Counter

PACKAGE = pathlib.Path(__file__).resolve().parents[1]
RESEARCH = PACKAGE.parent
QCODE = RESEARCH/'experimenting/qcode-open-entries'
BASE = QCODE/'general-criterion/unified'
sys.path[:0] = [str(BASE),str(BASE/'m8r1'),str(QCODE/'general')]

def sha(path):
    digest=hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):digest.update(block)
    return digest.hexdigest()

def read(name):return json.loads((BASE/name).read_text())

def check_hashes(obj):
    count=0
    for name,digest in obj.get('hashes',{}).items():
        p=pathlib.Path(name)
        if not p.is_absolute():p=BASE/p
        assert p.is_file() and sha(p)==digest, str(p)
        count+=1
    return count

def capture(module, filename, function='main'):
    """Redirect only the final report sink; all frozen data remain read-only."""
    m=importlib.import_module(module)
    target=PACKAGE/'audit'/filename
    saved=m.atomic_json
    def sink(path,value):
        target.write_text(json.dumps(value,indent=2)+'\n')
    m.atomic_json=sink
    try:
        with contextlib.redirect_stdout(io.StringIO()):getattr(m,function)()
    finally:m.atomic_json=saved
    return json.loads(target.read_text())

def ordinary_identity():
    """Independently check A6+A7<=2 using univariate MacWilliams rows."""
    from math import comb
    n=17;B=2**14
    def K(j,w):
        return sum((-1)**h*3**(j-h)*comb(w,h)*comb(n-w,j-h)
                   for h in range(max(0,j-(n-w)),min(j,w)+1))
    # Equality rows: A0=1, sum A=B, B_j=A_j for j<=5,
    # and A1=...=A5=0. Shadow rows are -signed-transform<=0.
    residual=[0]*18;rhs=0
    def add(mult,row,value=0):
        nonlocal rhs
        residual[:]=[a+mult*b for a,b in zip(residual,row)];rhs+=mult*value
    add(-1835008,[int(w==0) for w in range(18)],1)
    add(116,[1]*18,B)
    for j,mult in enumerate([60,25,11,3,1],1):
        add(mult,[K(j,w)-B*int(j==w) for w in range(18)])
    for j,mult in enumerate([-851968,-278528,-507904,-147456,-180224],1):
        add(mult,[int(w==j) for w in range(18)])
    for j,mult in enumerate([88,32,18,4,2]):
        add(mult,[-(-1)**w*K(j,w) for w in range(18)])
    assert residual==[32768*int(w in (6,7)) for w in range(18)]
    assert rhs==65536
    return {'coefficient':32768,'rhs':rhs,'conclusion':'A6+A7 <= 2'}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--deep',action='store_true')
    args=ap.parse_args();assert __debug__
    start=time.monotonic();report={'status':'RUNNING','checks':{}}
    audit=PACKAGE/'audit';audit.mkdir(exist_ok=True)
    def done(name,value):
        report['checks'][name]=value
        print(name, str(value)[:180],flush=True)
    certificate_files=sorted((QCODE/'pure/certificates_17_3_6').glob('*.json'))
    certificate_files += [BASE/'weight8_recenter/pure17_without_weight8.json']
    for group,checker,files in [
        ('16',QCODE/'release16/certificates/check_certificates.py',[]),
        ('17',QCODE/'pure/check_general.py',certificate_files)]:
        run=subprocess.run([sys.executable,'-B',str(checker),*map(str,files)],
                           text=True,capture_output=True,check=True)
        (audit/f'certificates-{group}.log').write_text(run.stdout)
        done('certificates_'+group,run.stdout.strip().splitlines()[-1])
    done('profile_identity',ordinary_identity())
    bundle=read('weight8_split/frame_count_certificates.json')
    assert len(bundle)==13
    with tempfile.TemporaryDirectory() as folder:
        paths=[]
        for name,data in bundle.items():
            p=pathlib.Path(folder)/name;p.write_text(json.dumps(data));paths.append(p)
        r=subprocess.run([sys.executable,'-B',str(QCODE/'pure/check_general.py'),*map(str,paths)],text=True,capture_output=True,check=True)
        (audit/'frame-count-certificates.log').write_text(r.stdout)
    done('frame_count_certificates',13)

    chain=capture('parent_mass_check.verify_chain','classification-chain.json')
    assert chain['status']=='PASS_EXHAUSTIVE_SUFFIX_CLASSIFICATION_CHAIN'
    done('classification_chain',chain['cells'])
    for name in ['prefix_mass_check/verification.json','kernel_case_plan/provenance.json',
                 'weight8_audit/audit.json','mixed_frames/stage3/profile_verification.json']:
        d=read(name);done(name,{'status':d['status'],'current_hashes':check_hashes(d)})

    # Full rank-two mixed frames: fresh read-only partition verification.
    from mixed_frames.stage2.common import hashes as stage2_hashes
    from mixed_frames.root_crosscheck.strict_aggregate import check_partition
    manifest=read('mixed_frames/stage2/manifest.json')
    assert stage2_hashes()==manifest['hashes']
    roots=leaves=0;classes=[]
    for row in manifest['classes']:
        if row['frame']!='668':continue
        for first in row['first_images']:
            folder=BASE/'mixed_frames/stage2/runs'/f'668_{row["class_index"]}'/str(first)
            cp=json.loads((folder/'checkpoint.json').read_text())
            for k,v in {'frame':'668','class_index':row['class_index'],'first_image':first,
                        'd':6,'mu_out':8,'hashes':manifest['hashes'],'input_sha256':row['input_sha256']}.items():assert cp[k]==v
            ok,why,survivors=check_partition(folder,row)
            assert ok and survivors==0,(folder,why,survivors)
            assert cp['excluded'] and cp['surviving_solutions']==0
            roots+=1;leaves+=cp['enumerator']['leaf_count']
        classes.append(row['class_index'])
    assert len(classes)==len(set(classes))==48
    done('frame_668',{'candidate_classes':48,'verified_roots':roots,'Z_leaves':leaves})
    for module,label in [('fast778.run','frame-778.json'),('speedup_review.run','frame-688.json')]:
        r=capture(module,label,'aggregate')
        assert r['status']=='COMPLETE' and r['verified'] and all(c['excluded'] for c in r['classes'])
        done(label,{'cases':len(r['classes']),'verified_roots':sum(c['partitions'] for c in r['classes'])})
    r=capture('rank2_kernel.combine_888','frame-888.json')
    assert r['status']=='CONDITIONALLY_COMPLETE_888_COLUMN' and not r['missing'] and r['classes']==56
    done('frame_888',{k:r[k] for k in ['classes','new_slices','new_Y_leaves_rejected','legacy_hull0_verified_roots']})

    for module,name,expected in [('radical_projection.aggregate','rankone-q3-hull5.json',118),
                                ('fixed_kernel.hull3_aggregate','rankone-q3-hull3.json',5995)]:
        r=capture(module,name)
        assert r['status']=='CONDITIONALLY_COMPLETE' and r['closed_classes']==r['classes']==expected
        assert all(not c.get('errors') for c in r['records'])
        done(name,{'closed_classes':expected,'counts':r['counts']})

    from oracle_residual.common import assignment,fingerprints
    from rankone_reuse.records import accepted_hashes
    from rankone_reuse.pool import validate
    from oracle_residual.aggregate import aggregate, verify as residual_verify
    from rankone_residual.search import expected_slices
    from rankone_reuse.records import record_path
    plan=validate();ass=assignment()
    saved=read('oracle_residual/aggregate.json')
    accepted=accepted_hashes()+[fingerprints()]
    assert saved['accepted_hashes']==accepted
    assert saved['aggregation_source_sha256']==sha(BASE/'oracle_residual/aggregate.py')
    assert saved['assignment_sha256']==sha(BASE/'oracle_residual/assignment.json')
    assert saved['status']=='CONDITIONALLY_COMPLETE_RESIDUAL_RANKONE'
    assert saved['selected_classes']==saved['classes']==len(saved['records'])==len(plan['records'])==32747
    assert not saved['errors'] and not saved['orphans']
    ids={r['cls'] for r in plan['records']}
    assert {r['cls'] for r in saved['records']}==ids
    for r in saved['records']:
        assert r['status']=='CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES'
        assert all(c['status']=='CONDITIONALLY_EXCLUDED_BRANCH' for c in r['branches'])
    assert saved['counts']['EXCLUDED_SLICE_CHUNK']==len(saved['record_hashes'])==2419796
    done('residual_saved_coverage',{'classes':32747,'records':2419796,'report_sha256':sha(BASE/'oracle_residual/aggregate.json'),'source_fingerprints_current':True})
    # Fresh sample includes both producer domains and every residual branch type.
    by_id={r['cls']:r for r in plan['records']}
    selected=set()
    for domain in [ass['local'],ass['remote']]:
        for hull in [1,3]:
            candidates=[i for i in domain if by_id[i]['hull']==hull]
            selected.update([min(candidates),max(candidates)])
    if args.deep:
        result=aggregate(plan,None)
        assert not result['errors'] and not result['orphans']
        assert all(c['status']=='CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES' for c in result['records'])
        sample={'status':result['status'],'selected_classes':len(result['records']),'counts':result['counts'],
                'record_hash_count':len(result['record_hashes']),'records':result['records']}
    else:
        checked=[]
        for cls in sorted(selected):
            row=by_id[cls];count=0
            for s in expected_slices(plan,row):
                for k,name in enumerate(s['names']):
                    if cls in set(ass['remote']):
                        path=BASE/'oracle_residual/runs'/name
                    else:
                        path=record_path(name)
                    assert path is not None and path.is_file(),name
                    assert residual_verify(path,row,s,k,accepted=accepted)=='EXCLUDED_SLICE_CHUNK'
                    count+=1
            checked.append({'cls':cls,'records_checked':count})
            print('fresh_residual_class',cls,count,flush=True)
        sample={'status':'PASS_COMPLETE_SELECTED_CLASS_REPLAYS','selected_classes':len(checked),
                'record_hash_count':sum(r['records_checked'] for r in checked),'records':checked}
    (audit/'residual-fresh-check.json').write_text(json.dumps(sample,indent=2)+'\n')
    done('residual_fresh_check',{'deep':args.deep,'classes':sample['selected_classes'],'records':sample['record_hash_count']})
    census=read('opus_parallel/suffix9_dimension5/suffix9_dim5.json')['classes']
    hulls=Counter(r['representative']['hull'] for r in census)
    assert hulls==Counter({1:26752,3:5995,5:118})
    done('rank_one_case_exhaustiveness',{'suffix_classes':32865,'hulls':dict(hulls),'kernel_branches':[[1,0,1],[2,2,1],[2,0,3],[3,2,3],[3,0,5]]})
    report.update(status='PASS',elapsed_seconds=time.monotonic()-start,
        scope='Exact certificate replay, source-provenance validation, exhaustive saved nonresidual partition reaggregation, all residual class coverage and '+('full residual replay.' if args.deep else 'fresh eight-class residual replay. The full residual aggregation is preserved and hash-bound, not independently rerun in this default check.'))
    (audit/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('UNIFIED EVIDENCE AUDIT PASS',flush=True)

if __name__=='__main__':main()
