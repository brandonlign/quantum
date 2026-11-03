"""Complete heavy-plane partition of target17 m8,r1,hull5 extensions.

This branch forces q=3,t=0,dim(S intersect W)=2. Fixing that heavy plane
restricts all syndrome images to its annihilator, with only 2048 labels.
LIMIT in any plane prevents a complete class exclusion.
"""
import argparse
import hashlib
import json
import sys
import time
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from model import ExtensionModel, coordinates, gram, lift_coordinate, reduced_basis, weight
from run_benchmark import z_condition
from signed_extension import solve_extension
from run_m8r1 import leaves


def heavy_planes(C):
    words=[lift_coordinate(v,C) for v in range(32)]
    heavy={i for i,w in enumerate(words) if i and weight(w,9)>=8}
    return [(words[a],words[b]) for a in sorted(heavy) for b in sorted(heavy)
            if a<b<(a^b) and (a^b) in heavy]


def quotient_input(model,N):
    Ncoords=[coordinates(n,model.P) for n in N]
    assert len(reduced_basis(Ncoords))==2
    full=list(Ncoords);quotient=[]
    for i in range(model.b):
        v=1<<i
        if len(reduced_basis(full+[v]))>len(full):
            full.append(v);quotient.append(v)
    assert len(full)==13 and len(quotient)==11
    # Evaluations on [N,Pquotient] have first two values zero. Convert them
    # back to evaluations on the original ordered physical P basis.
    units=[sum(((coordinates(1<<i,full)>>(2+j))&1)<<i for i in range(13))
           for j in range(11)]
    labels=[lift_coordinate(v,units) for v in range(2048)]
    Q=[lift_coordinate(v,model.P) for v in quotient]
    return dict(N=list(N),Ncoords=Ncoords,units=units,labels=labels,
                distances=[model.distances[v] for v in labels],Gamma=gram(Q,9))


def run_plane(model,N,seconds,solver_nodes,backend='rank_xdom',first_representatives=None):
    data=quotient_input(model,N);stats={};tally=Counter();found=None;limits=0
    orbit_data=None
    if backend=='symmetry_h5':
        from suffix_symmetries import syndrome_orbits
        orbit_data=syndrome_orbits(model,data,N)
    if first_representatives is not None:
        if orbit_data is None or len(first_representatives)!=len(set(first_representatives)):
            raise ValueError('Partitions require symmetry backend and distinct verified representatives')
        if not set(first_representatives)<=set(orbit_data['first_representatives']):
            raise ValueError('Unverified first-image partition')
    selected=None if orbit_data is None else (
        orbit_data['first_representatives'] if first_representatives is None else first_representatives)
    for small in leaves(data['distances'],seconds,None,stats,
                        Gamma=data['Gamma'],backend=backend,
                        first_representatives=selected):
        xs=[lift_coordinate(v,data['units']) for v in small]
        assert z_condition(xs,model.Gamma,13)
        answer=solve_extension(model,xs,node_limit=solver_nodes,minimum_stabilizer_weight=8)
        tally[answer['status']+':'+answer.get('reason','')]+=1
        if answer['status']=='LIMIT':limits+=1
        if answer['status']=='FOUND':found=dict(xs=xs,answer=answer);break
    status='FOUND' if found else ('EXCLUDED' if stats.get('status')=='COMPLETE' and not limits else 'LIMIT')
    if status=='EXCLUDED' and first_representatives is not None:
        status='PARTITION_EXCLUDED'
    return dict(status=status,N=list(N),z=stats,tally=dict(tally),found=found,orbit_data=orbit_data,
                searched_first_representatives=selected)


def run_class(C,seconds,solver_nodes,backend='rank_xdom'):
    start=time.monotonic();model=ExtensionModel(17,3,6,8,[255<<8],C)
    assert model.hull==5 and model.C_distance>=6
    planes=heavy_planes(model.C)
    results=[]
    for plane in planes:
        r=run_plane(model,plane,seconds,solver_nodes,backend);results.append(r)
        if r['status']=='FOUND':break
    status='FOUND' if any(r['status']=='FOUND' for r in results) else (
        'EXCLUDED' if all(r['status']=='EXCLUDED' for r in results) else 'LIMIT')
    return dict(status=status,backend=backend,plane_count=len(planes),planes=results,
                seconds=round(time.monotonic()-start,3))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--only',required=True)
    ap.add_argument('--seconds',type=float,default=15)
    ap.add_argument('--solver-nodes',type=int,default=10000)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--backend',choices=['rank_xdom','fixed_h5','avoid_h5','symmetry_h5'],default='rank_xdom')
    a=ap.parse_args()
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    classes=json.loads(source.read_text())['classes'];records=[]
    from run_m8r1 import BACKENDS
    binary,_=BACKENDS[a.backend]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [source,HERE/'run_hull5_planes.py',HERE/'run_m8r1.py',HERE/binary,HERE/(binary+'.c')]}
    if a.backend=='symmetry_h5':
        p=HERE/'suffix_symmetries.py';hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    for index in map(int,a.only.split(',')):
        result=dict(cls=index,**run_class(classes[index]['representative']['generators'],a.seconds,a.solver_nodes,a.backend))
        records.append(result)
        print(json.dumps(result),flush=True)
        a.output.write_text(json.dumps(dict(records=records,sha256=hashes,
            scope='Hull-five class only; exhausted heavy-plane partition and all fixed-Z problems required'),indent=2)+'\n')
