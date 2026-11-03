"""Combine complete class exclusions, preserving provenance and scope."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from run_hull5_planes import heavy_planes,reduced_basis
from heavy_plane_orbits import geometry

HERE=Path(__file__).resolve().parent


def main():
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    classes=json.loads(source.read_text())['classes']
    ids={c['class_index'] for c in classes if c['representative']['hull']==5}
    proof={};inputs=[source]
    def add(index,why):
        assert index in ids
        proof.setdefault(index,[]).append(why)
    for name in ['local_hull5_plane_benchmark.json','single_plane_batch.json','single_plane_retry.json',
                 'single_plane_retry_remaining.json','class29895_symmetry.json']:
        p=HERE/name
        if not p.exists():continue
        inputs.append(p)
        for r in json.loads(p.read_text())['records']:
            if r['status']=='EXCLUDED':add(r['cls'],name)
    p=HERE/'class29895_partitions.json'
    if p.exists():
        d=json.loads(p.read_text());inputs.append(p)
        if d['status']=='EXCLUDED':
            v=json.loads((HERE/'partition_coverage_verification.json').read_text())
            assert v['status']=='PASS' and v['sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
            add(d['cls'],'verified complete first-orbit partitions')
            inputs.append(HERE/'partition_coverage_verification.json')
    p=HERE.parent/'projection_rank_verification.json';d=json.loads(p.read_text());inputs.append(p)
    assert d['status']=='PASS'
    for r in d['target_analytic_exclusions']:add(r['cls'],'general projection-rank obstruction')
    r2path=HERE.parent/'rank2_branch_result.json';r2=json.loads(r2path.read_text());inputs.append(r2path)
    assert r2['status']=='COMPLETE_BRANCH_EXCLUSION'
    assert r2['local_branch']=={'m':8,'r':2}
    for c in classes:
        if c['class_index'] not in ids:continue
        pp=heavy_planes(reduced_basis(c['representative']['generators']))
        if pp and all(geometry(*p)=='common_zero' for p in pp):
            add(c['class_index'],'completed m8,r2 case propagated to every heavy plane')
    p=HERE/'full_support_geometry_audit.json';audit=json.loads(p.read_text());inputs.append(p)
    for r in audit['records']:
        if r['status']=='EXCLUDED_ALTERNATE':add(r['cls'],'alternate full-support-geometry search audit')
    # Different checkpoint files can prove the same class; count once.
    p=HERE/'multi_plane_pilot.json'
    if p.exists():
        d=json.loads(p.read_text());inputs.append(p)
        for index,status in d['statuses'].items():
            if status=='EXCLUDED':add(int(index),'complete verified heavy-plane orbit coverage')
        # Separately verified resumable first-image partitions can finish
        # individual plane orbits; never treat one plane as a whole class.
        closed_orbits={key for key,r in d['records'].items() if r['status'] in {
            'EXCLUDED','EXCLUDED_RANK_TWO_PROPAGATION','EXCLUDED_FULL_SUPPORT_AUDIT'}}
        for verified_path in HERE.glob('class*_plane*_partitions_verified.json'):
            v=json.loads(verified_path.read_text())
            checkpoint=verified_path.with_name(verified_path.name.replace('_verified.json','.json'))
            assert v['status']=='PASS' and v['sha256']==hashlib.sha256(checkpoint.read_bytes()).hexdigest()
            partition=json.loads(checkpoint.read_text())
            assert partition['status']=='EXCLUDED'
            for name,digest in partition['sha256'].items():
                assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest
            closed_orbits.add(f"{v['cls']}:{v['plane_orbit']}")
            inputs.extend([verified_path,checkpoint])
        plane_plan=HERE/'heavy_plane_orbits.json';inputs.append(plane_plan)
        for c in json.loads(plane_plan.read_text())['records']:
            if all(f"{c['cls']}:{i}" in closed_orbits for i in range(len(c['orbits']))):
                add(c['cls'],'complete plane coverage including independently audited first-image partitions')
    singles={c['class_index'] for c in classes if c['class_index'] in ids and
             len(heavy_planes(reduced_basis(c['representative']['generators'])))==1}
    result=dict(status='PARTIAL_HULL5' if len(proof)<len(ids) else 'COMPLETE_HULL5_EXCLUSION',
                hull5_classes=len(ids),excluded_classes=len(proof),remaining_classes=sorted(ids-set(proof)),
                single_plane_classes=len(singles),excluded_single_plane_classes=len(singles.intersection(proof)),
                single_plane_subcase_complete=singles<=set(proof),
                provenance={str(i):why for i,why in sorted(proof.items())},
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                scope='Hull-five portion of m8,r1 only; other hulls and both other rank-one branches remain unresolved')
    (HERE/'hull5_aggregate.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['provenance','sha256','remaining_classes']})


if __name__=='__main__':main()
