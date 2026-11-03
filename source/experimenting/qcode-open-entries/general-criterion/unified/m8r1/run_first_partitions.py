"""Resume symmetry-first-image partitions without redoing completed pieces."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
from run_hull5_planes import ExtensionModel,heavy_planes,quotient_input,run_plane
from suffix_symmetries import syndrome_orbits

HERE=Path(__file__).resolve().parent


def aggregate(records,required):
    covered=set()
    if any(r['status']=='FOUND' for r in records):return 'FOUND',sorted(covered)
    for r in records:
        if r['status']=='PARTITION_EXCLUDED':covered.update(r['searched_first_representatives'])
    return ('EXCLUDED' if covered==set(required) else 'PARTIAL'),sorted(covered)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--class',dest='index',type=int,default=29895)
    ap.add_argument('--seconds',type=float,default=30);ap.add_argument('--solver-nodes',type=int,default=10000)
    ap.add_argument('--workers',type=int,default=2);ap.add_argument('--chunk',type=int,default=10)
    ap.add_argument('--output',type=Path,default=HERE/'class29895_partitions.json')
    a=ap.parse_args()
    if a.workers not in [1,2] or a.chunk<1:raise ValueError('Use1/2 local CPU jobs and positive chunk size')
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    classes=json.loads(source.read_text())['classes']
    model=ExtensionModel(17,3,6,8,[255<<8],classes[a.index]['representative']['generators'])
    assert model.hull==5
    planes=heavy_planes(model.C)
    if len(planes)!=1:raise ValueError('This class-level aggregator requires exactly one heavy plane')
    N=planes[0];data=quotient_input(model,N);orbits=syndrome_orbits(model,data,N)
    reps=orbits['first_representatives']
    # Required subgroup-orbit coverage is fixed and checked on every resume.
    paths=[source,HERE/'z_enum_symmetry.c',HERE/'z_enum_symmetry',HERE/'suffix_symmetries.py',
           HERE/'run_hull5_planes.py',HERE/'run_m8r1.py',Path(__file__)]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    records=[]
    if a.output.exists():
        old=json.loads(a.output.read_text())
        if old['sha256']!=hashes or old['cls']!=a.index or old['required_first_representatives']!=reps:
            raise ValueError('Resume inputs/code changed; audit before reusing completed partitions')
        records=old['records']
    status,covered=aggregate(records,reps)
    remaining=[r for r in reps if r not in covered]
    groups=[remaining[i:i+a.chunk] for i in range(0,len(remaining),a.chunk)]
    def save():
        status,covered=aggregate(records,reps)
        result=dict(cls=a.index,status=status,N=list(N),required_first_representatives=reps,
                    excluded_first_representatives=covered,records=records,sha256=hashes,
                    scope='One heavy plane; class exclusion requires exhaustive coverage of every verified first-image orbit and resolved fixed-Z searches')
        tmp=a.output.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(a.output)
    save()
    if status in ['EXCLUDED','FOUND']:print(status,flush=True);return
    print(json.dumps(dict(cls=a.index,remaining_first=len(remaining),groups=len(groups))),flush=True)
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        jobs={pool.submit(run_plane,model,N,a.seconds,a.solver_nodes,'symmetry_h5',group):group for group in groups}
        for job in as_completed(jobs):
            record=job.result();records.append(record);save()
            status,covered=aggregate(records,reps)
            print(json.dumps(dict(first_group=jobs[job],status=record['status'],
                                  nodes=record['z'].get('nodes'),leaves=record['z'].get('leaves'),
                                  covered=len(covered),total=len(reps))),flush=True)
    print(json.dumps(dict(status=aggregate(records,reps)[0],covered=len(aggregate(records,reps)[1]))),flush=True)


if __name__=='__main__':main()
