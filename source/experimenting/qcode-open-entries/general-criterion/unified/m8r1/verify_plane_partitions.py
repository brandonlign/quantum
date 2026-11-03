"""Verify complete first-orbit coverage, status semantics and input hashes."""
import argparse
import hashlib
import json
from pathlib import Path
from run_hull5_planes import ExtensionModel,heavy_planes,quotient_input
from suffix_symmetries import syndrome_orbits
from run_plane_first_partitions import aggregate

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('path',type=Path);a=ap.parse_args()
    path=a.path;d=json.loads(path.read_text())
    assert d['status']=='EXCLUDED'
    for name,digest in d['sha256'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    C=json.loads(source.read_text())['classes'][d['cls']]['representative']['generators']
    model=ExtensionModel(17,3,6,8,[255<<8],C);planes=heavy_planes(model.C)
    assert tuple(d['N']) in planes
    plan=json.loads((HERE/'heavy_plane_orbits.json').read_text())['records']
    selected=next(r for r in plan if r['cls']==d['cls'])['orbits'][d['plane_orbit']]
    assert selected['representative']==d['N']
    orbits=syndrome_orbits(model,quotient_input(model,d['N']),d['N'])
    required=orbits['first_representatives']
    assert required==d['required_first_representatives']
    covered=set();complete=0;nodes=leaves=0
    for r in d['records']:
        if r['status']!='PARTITION_EXCLUDED':continue
        assert r['z']['status']=='COMPLETE' and r['found'] is None
        assert all(key.startswith('EXCLUDED:') for key in r['tally'])
        assert sum(r['tally'].values())==r['z']['leaves']
        first=r['searched_first_representatives']
        assert set(first)<=set(required) and not covered.intersection(first)
        covered.update(first);complete+=1;nodes+=r['z']['nodes'];leaves+=r['z']['leaves']
    assert covered==set(required) and sorted(covered)==d['excluded_first_representatives']
    # Missing, duplicate, or solver-limited partitions cannot close a class.
    def rec(first,status='PARTITION_EXCLUDED'):return dict(status=status,searched_first_representatives=first)
    assert aggregate([rec([3])],[3,12])[0]=='PARTIAL'
    assert aggregate([rec([3]),rec([3])],[3,12])[0]=='PARTIAL'
    assert aggregate([rec([3]),rec([12],'LIMIT')],[3,12])[0]=='PARTIAL'
    assert aggregate([rec([3]),rec([12])],[3,12])[0]=='EXCLUDED'
    result=dict(status='PASS',cls=d['cls'],plane_orbit=d['plane_orbit'],required_first_orbits=len(required),complete_partitions=complete,
                proof_nodes=nodes,proof_leaves=leaves,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Coverage/input/status audit with physically rederived orbits; no second exhaustive rerun of all search nodes')
    path.with_name(path.stem+'_verified.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)


if __name__=='__main__':main()
