"""Rebuild first-image orbits and audit complete partition coverage/statuses."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from model import lift_coordinate
from rank2_mu7.run_pair_partitions import aggregate

BASE=Path(__file__).resolve().parent.parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('checkpoint',type=Path);a=ap.parse_args()
    d=json.loads(a.checkpoint.read_text());assert d['status']=='EXCLUDED_MIN7_CANDIDATE'
    for name,digest in d['sha256'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest
    index=d['class_index'];table_path=BASE/f'rank2_z/class_{index}.json'
    table=json.loads(table_path.read_text())
    auto=next(r for r in json.loads((BASE/'rank2_z/aut_orbits.json').read_text())['classes'] if r['class_index']==index)
    assert auto['table_sha256']==hashlib.sha256(table_path.read_bytes()).hexdigest()
    high={i for i,(delta,lam) in enumerate(zip(table['coset_distances'],table['suffix_weights'])) if i and delta>=4 and lam>=5}
    remaining=set(high);required=[]
    while remaining:
        seed=min(remaining);orbit={seed};front=[seed]
        while front:
            v=front.pop()
            for T in auto['verified_generator_coefficient_maps']:
                w=lift_coordinate(v,T);assert w in high
                if w not in orbit:orbit.add(w);front.append(w)
        assert orbit<=remaining;remaining-=orbit;required.append(seed)
    assert d['required_first']==required
    covered=set();nodes=leaves=0
    for r in d['records']:
        if r['status']!='PARTITION_EXCLUDED':continue
        assert r['z']['status']=='COMPLETE' and r['found'] is None
        assert r['z']['leaves']==len(r['verdicts'])
        assert all(v['answer']['status']=='EXCLUDED' for v in r['verdicts'])
        assert len(r['first'])==len(set(r['first'])) and set(r['first'])<=set(required)
        assert not covered.intersection(r['first']);covered.update(r['first'])
        path=Path(r['input']);log=Path(r['log'])
        assert hashlib.sha256(path.read_bytes()).hexdigest()==r['input_sha256']
        assert hashlib.sha256(log.read_bytes()).hexdigest()==r['log_sha256']
        rows=path.read_text().splitlines();assert list(map(int,rows[-1].split()))==r['first']
        assert int(rows[-2])==len(r['first']) and rows[0].split()[-1]=='0'
        base_rows=(BASE/f'rank2_mu7_orbits/class_{index}.txt').read_text().splitlines()
        assert rows[1:4098]==base_rows[1:4098]
        events=[json.loads(line) for line in log.read_text().splitlines()]
        assert events[-1]==r['z']
        assert [e['Z'] for e in events if e['status']=='LEAF']==[v['Z'] for v in r['verdicts']]
        nodes+=r['z']['nodes'];leaves+=r['z']['leaves']
    assert sorted(covered)==required==d['covered_first']
    def rec(first,status='PARTITION_EXCLUDED'):return dict(status=status,first=first)
    assert aggregate([rec([3])],[3,5])[0]=='PARTIAL'
    assert aggregate([rec([3]),rec([3])],[3,5])[0]=='PARTIAL'
    assert aggregate([rec([3]),rec([5],'LIMIT')],[3,5])[0]=='PARTIAL'
    assert aggregate([rec([3]),rec([5])],[3,5])[0]=='EXCLUDED_MIN7_CANDIDATE'
    result=dict(status='PASS',class_index=index,required_first=len(required),
                complete_partitions=sum(r['status']=='PARTITION_EXCLUDED' for r in d['records']),
                proof_nodes=nodes,proof_leaves=leaves,checkpoint_sha256=hashlib.sha256(a.checkpoint.read_bytes()).hexdigest(),
                scope='Orbit rederivation, disjoint complete coverage, inputs/logs/statuses checked; no exhaustive node replay or independent automorphism implementation')
    a.checkpoint.with_name(a.checkpoint.stem+'_verified.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)


if __name__=='__main__':main()
