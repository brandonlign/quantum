"""Resume hull-five plane orbits, with explicit cross-case dependencies."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
from run_hull5_planes import ExtensionModel,run_plane
from heavy_plane_orbits import geometry

HERE=Path(__file__).resolve().parent
FINAL={'EXCLUDED','EXCLUDED_RANK_TWO_PROPAGATION','EXCLUDED_FULL_SUPPORT_AUDIT'}


def statuses(plan,records):
    out={}
    for c in plan:
        keys=[f"{c['cls']}:{i}" for i in range(len(c['orbits']))]
        values=[records.get(k,{}).get('status') for k in keys]
        out[str(c['cls'])]='FOUND' if 'FOUND' in values else (
            'EXCLUDED' if all(v in FINAL for v in values) else 'PARTIAL')
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--only',required=True)
    ap.add_argument('--seconds',type=float,default=30);ap.add_argument('--solver-nodes',type=int,default=10000)
    ap.add_argument('--workers',type=int,default=2);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.workers not in [1,2]:raise ValueError('Use1/2 local CPU jobs')
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    orbit_path=HERE/'heavy_plane_orbits.json'
    r2path=HERE.parent/'rank2_branch_result.json';audit_path=HERE/'full_support_geometry_audit.json'
    classes=json.loads(source.read_text())['classes'];plan=json.loads(orbit_path.read_text())['records']
    r2=json.loads(r2path.read_text());audit=json.loads(audit_path.read_text())
    assert r2['status']=='COMPLETE_BRANCH_EXCLUSION' and r2['local_branch']=={'m':8,'r':2}
    audited={r['cls'] for r in audit['records'] if r['status']=='EXCLUDED_ALTERNATE'}
    wanted=set(map(int,a.only.split(',')));selected=[r for r in plan if r['cls'] in wanted]
    assert len(selected)==len(wanted)
    paths=[source,orbit_path,r2path,audit_path,Path(__file__),HERE/'run_hull5_planes.py',
           HERE/'run_m8r1.py',HERE/'suffix_symmetries.py',HERE/'z_enum_symmetry.c',HERE/'z_enum_symmetry']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    records={}
    if a.output.exists():
        old=json.loads(a.output.read_text())
        if old['sha256']!=hashes:raise ValueError('Resume sources/code changed; audit before reusing records')
        records=old['records']
    def save():
        data=dict(statuses=statuses(plan,records),records=records,sha256=hashes,
                  scope='Class exclusion requires all verified heavy-plane orbits closed; LIMIT/missing orbit is unresolved')
        tmp=a.output.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(a.output)
    jobs=[];models={}
    for c in selected:
        model=ExtensionModel(17,3,6,8,[255<<8],classes[c['cls']]['representative']['generators'])
        assert model.hull==5;models[c['cls']]=model
        for i,o in enumerate(c['orbits']):
            key=f"{c['cls']}:{i}"
            if records.get(key,{}).get('status') in FINAL:continue
            N=o['representative'];kind=geometry(*N)
            assert kind==o['geometry']
            if kind=='common_zero':
                records[key]=dict(status='EXCLUDED_RANK_TWO_PROPAGATION',N=N,orbit_size=o['size'])
            elif kind=='full_support_word':
                assert c['cls'] in audited
                records[key]=dict(status='EXCLUDED_FULL_SUPPORT_AUDIT',N=N,orbit_size=o['size'])
            else:jobs.append((key,c['cls'],N,o['size']))
    save()
    print(json.dumps(dict(selected_classes=len(selected),pending_plane_orbits=len(jobs))),flush=True)
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        pending={pool.submit(run_plane,models[index],N,a.seconds,a.solver_nodes,'symmetry_h5'):
                 (key,index,N,size) for key,index,N,size in jobs}
        for job in as_completed(pending):
            key,index,N,size=pending[job];r=job.result();r['orbit_size']=size;records[key]=r;save()
            print(json.dumps(dict(key=key,status=r['status'],nodes=r['z'].get('nodes'),leaves=r['z'].get('leaves'),
                                  class_status=statuses(plan,records)[str(index)])),flush=True)
    print(json.dumps({str(c['cls']):statuses(plan,records)[str(c['cls'])] for c in selected}),flush=True)


if __name__=='__main__':main()
