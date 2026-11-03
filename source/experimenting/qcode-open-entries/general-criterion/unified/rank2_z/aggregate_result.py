"""Consolidate only completed, input-hash-verified class exclusions."""
from pathlib import Path
import json,hashlib
BASE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    inputs=json.loads((BASE.parent/'rank2_z_inputs.json').read_text());expected={v['class_index'] for v in inputs['classes']}
    tables=json.loads((BASE/'table_verification.json').read_text());autos=json.loads((BASE/'aut_orbits.json').read_text())
    assert tables['status']==autos['status']=='PASS'
    backend=json.loads((BASE/'all_forward_verification.json').read_text());assert backend['status']=='PASS'
    records=[];sources=[BASE/'symmetric_long.json',BASE/'remaining_long.json']
    for source in sources:
        batch=json.loads(source.read_text());assert batch['symmetry_first_orbits']
        assert batch['sha256']['enumerator_all_unary.cpp']==sha(BASE/'enumerator_all_unary.cpp')
        assert batch['sha256']['enumerator_all_unary']==sha(BASE/'enumerator_all_unary')
        for r in batch['reports']:
            assert r['status']=='COMPLETE' and r['Z_excluded'] and r['leaves']==0
            log=Path(r['log']);end=json.loads(log.read_text().splitlines()[-1])
            assert end['status']=='COMPLETE' and end['Z_excluded'] and end['leaves']==0
            txt=log.with_suffix('.txt');assert sha(txt)==r['input_sha256']
            record=dict(r,report_source=str(source),log_sha256=sha(log))
            records.append(record)
    assert len(records)==len(expected)==15 and {r['class_index'] for r in records}==expected
    for item in inputs['classes']:
        matching=[v for v in tables['reports'] if v['class_index']==item['class_index']][0]
        matching_auto=[v for v in autos['classes'] if v['class_index']==item['class_index']][0]
        assert matching['table_sha256']==matching_auto['table_sha256']==sha(Path(item['table']))
    suffix=json.loads((BASE.parent/'suffix9_verification.json').read_text())
    equivalence=json.loads((BASE.parent/'suffix9_equivalence_verification.json').read_text())
    assert suffix['status']==equivalence['status']=='PASS'
    result=dict(status='COMPLETE_BRANCH_EXCLUSION',target='[[17,3,6]]',local_branch=dict(m=8,r=2),
        initial_suffix_classes=56,nondegenerate_candidates=15,excluded_candidates=15,
        total_search_nodes=sum(r['nodes'] for r in records),
        search_seconds=sum(r['seconds'] for r in records),records=sorted(records,key=lambda r:r['class_index']),
        scope='Complete exclusion of the m8,r2 local branch, conditional on the supplied '
              'mass-checked parent census and the documented exact extension reduction. '
              'Three rank-one local branches remain; the full parameter set is unresolved.',
        proof_note=str(BASE.parent/'RANK_TWO_PROGRESS.md'),
        sha256={str(p):sha(p) for p in [Path(__file__),*sources,BASE/'table_verification.json',
            BASE/'aut_orbits.json',BASE/'all_forward_verification.json',
            BASE.parent/'suffix9_verification.json',BASE.parent/'suffix9_equivalence_verification.json']})
    (BASE.parent/'rank2_branch_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['records','sha256']}))
if __name__=='__main__':main()
