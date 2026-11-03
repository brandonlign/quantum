"""Alternate full-Gram/radical exhaustive search and direct leaf rejection.

The distance-bitset/canonical-enumeration core is shared with the first
backend. Commutation/radical pruning uses a different calculation, and
unary-X forward pruning and signed_extension are not used.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from run_hull5_planes import quotient_input,heavy_planes,ExtensionModel,leaves
from verify_hull5_rejections import direct_rejection
from model import lift_coordinate

HERE=Path(__file__).resolve().parent

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--only',default='150,181,217')
    ap.add_argument('--seconds',type=float,default=90)
    ap.add_argument('--output',type=Path,default=HERE/'hull5_alternate_verification.json')
    a=ap.parse_args()
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    classes=json.loads(source.read_text())['classes'];records=[]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [source,HERE/'z_enum_avoid_radical.c',HERE/'z_enum_avoid_radical',Path(__file__)]}
    for index in map(int,a.only.split(',')):
        t0=time.monotonic();model=ExtensionModel(17,3,6,8,[255<<8],classes[index]['representative']['generators'])
        planes=[]
        for N in heavy_planes(model.C):
            data=quotient_input(model,N);stats={};rejected=[]
            for small in leaves(data['distances'],a.seconds,None,stats,Gamma=data['Gamma'],backend='avoid_h5'):
                xs=[lift_coordinate(v,data['units']) for v in small]
                rejected.append(dict(xs=xs,**direct_rejection(model,xs)))
            planes.append(dict(N=list(N),stats=stats,rejections=rejected))
        result=dict(cls=index,status='EXCLUDED_ALTERNATE' if all(p['stats']['status']=='COMPLETE' for p in planes) else 'LIMIT',
                    planes=planes,seconds=round(time.monotonic()-t0,3))
        records.append(result)
        print(json.dumps({k:v for k,v in result.items() if k!='planes'}),flush=True)
        a.output.write_text(json.dumps(dict(records=records,sha256=hashes,
            scope='Alternate full-annihilator/radical pruning with direct physical necessary leaf rejection; canonical distance bitset core shared'),indent=2)+'\n')
