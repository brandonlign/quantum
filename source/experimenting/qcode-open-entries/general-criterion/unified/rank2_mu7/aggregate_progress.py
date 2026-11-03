"""Combine mu7 exclusions without borrowing any minimum-eight result."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=HERE.parent


def main():
    candidates={r['class_index'] for r in json.loads((BASE/'rank2_z_inputs.json').read_text())['classes']}
    assert len(candidates)==15
    proof={};inputs=[]
    def add(index,why):
        assert index in candidates;proof.setdefault(index,[]).append(why)
    for name in ['rank2_mu7_orbits','rank2_mu7_pairs']:
        path=BASE/name/'screen.json'
        if not path.exists():continue
        data=json.loads(path.read_text());assert data['minimum_stabilizer']==7
        for file,digest in data['sha256'].items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest
        inputs.append(path)
        for r in data['records']:
            log=Path(r['log']);table=BASE/f"rank2_z/class_{r['class_index']}.json"
            path_input=log.with_suffix('.txt')
            assert hashlib.sha256(path_input.read_bytes()).hexdigest()==r['input_sha256']
            assert hashlib.sha256(table.read_bytes()).hexdigest()==r['source_sha256']
            events=[json.loads(line) for line in log.read_text().splitlines()]
            assert events[-1]['status']==r['status'] and events[-1]['leaves']==r['leaves']
            if r['status']=='COMPLETE' and r['leaves']==0:add(r['class_index'],name+' complete zero-leaf screen')
    for verified in HERE.glob('class*_verified.json'):
        audit=json.loads(verified.read_text());checkpoint=verified.with_name(verified.name.replace('_verified.json','.json'))
        data=json.loads(checkpoint.read_text())
        assert audit['status']=='PASS' and audit['checkpoint_sha256']==hashlib.sha256(checkpoint.read_bytes()).hexdigest()
        assert data['status']=='EXCLUDED_MIN7_CANDIDATE' and data['minimum_stabilizer']==7
        for file,digest in data['sha256'].items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest
        add(data['class_index'],'verified complete first-image partitions')
        inputs.extend([verified,checkpoint])
    result=dict(status='PARTIAL_RECENTERED_MIN7' if set(proof)!=candidates else 'COMPLETE_RECENTERED_MIN7_EXCLUSION',
                candidate_count=15,excluded_candidates=len(proof),remaining_candidates=sorted(candidates-set(proof)),
                provenance={str(k):v for k,v in sorted(proof.items())},
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                scope='Only pure[[17,3,6]] with minimum stabilizer>=7 and a full rank-two eight-site plane; other target branches remain separate')
    (HERE/'aggregate.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['provenance','sha256']})


if __name__=='__main__':main()
