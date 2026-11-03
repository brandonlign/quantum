"""Replay completed class searches; reject leaves with direct necessary tests.

Does not call signed_extension.solve_extension or its normal-form code.
The Z backend itself is shared; this is not a second independent full census.
"""
import hashlib
import json
from pathlib import Path
from run_hull5_planes import quotient_input,heavy_planes,ExtensionModel,leaves
from model import nullspace,lift_coordinate,span,weight

HERE=Path(__file__).resolve().parent


def direct_rejection(model,xs):
    # For any compatible map, the suffix-only stabilizer is exactly
    # rad(Gamma) intersect Ann(xs), since the xs are independent.
    forced=nullspace(model.Gamma+list(xs),13)
    words=[lift_coordinate(c,model.P) for c in span(forced)]
    bad=next((w for w in words if w and weight(w,9)<8),None)
    if bad is not None:
        return dict(reason='forced_suffix_stabilizer_weight',word=bad,weight=weight(bad,9))
    # Choose alpha_i dual to the Z columns. The matrix identity implies
    # y_i = Gamma alpha_i + span(xs), independently of all other y_j.
    Z=span(xs)
    translations=[lift_coordinate(v,xs) for v in range(128)]
    for i in range(7):
        alpha=next(a for a in range(8192)
                   if sum(((a&x).bit_count()%2)<<j for j,x in enumerate(xs))==1<<i)
        y0=sum(((row&alpha).bit_count()%2)<<j for j,row in enumerate(model.Gamma))
        px=(1<<i)|(1<<7)
        costs=[min(weight(px|(z<<8),8),weight(px|((z^255)<<8),8)) for z in range(128)]
        witnesses=[]
        for zshift in Z:
            y=y0^zshift
            failure=next((z for z in range(128)
                          if costs[z]+model.distances[y^translations[z]]<6),None)
            if failure is None:break
            witnesses.append(dict(y=y,z=failure,cost=costs[failure],
                                  delta=model.distances[y^translations[failure]]))
        else:
            assert len(witnesses)==128
            return dict(reason='empty_full_affine_unary_domain',X_site=i,
                        affine_offset=y0,witnesses=witnesses)
    raise AssertionError('Leaf needs a stronger independent rejection')


def verify():
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    classes=json.loads(source.read_text())['classes'];records=[]
    saved=json.loads((HERE/'local_hull5_plane_benchmark.json').read_text())
    for earlier in saved['records']:
        index=earlier['cls'];model=ExtensionModel(17,3,6,8,[255<<8],classes[index]['representative']['generators'])
        plane_records=[]
        for N in heavy_planes(model.C):
            data=quotient_input(model,N);stats={};rejects=[]
            for small in leaves(data['distances'],45,None,stats,Gamma=data['Gamma'],backend='rank_xdom'):
                xs=[lift_coordinate(v,data['units']) for v in small]
                rejects.append(dict(xs=xs,**direct_rejection(model,xs)))
            assert stats['status']=='COMPLETE'
            plane_records.append(dict(N=list(N),stats=stats,rejections=rejects))
        assert len(plane_records)==earlier['plane_count']
        assert sum(p['stats']['leaves'] for p in plane_records)==sum(p['z']['leaves'] for p in earlier['planes'])
        records.append(dict(cls=index,planes=plane_records,status='REPLAYED_EXCLUSION'))
    result=dict(status='PASS',records=records,
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in [source,HERE/'z_enum_rank.c',HERE/'z_enum_rank',Path(__file__)]},
                scope='Complete rerun with shared Z backend; independent direct stabilizer/unary leaf rejections. No independent second full Z enumeration.')
    (HERE/'hull5_rejection_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',classes=[r['cls'] for r in records],
                         leaves=sum(p['stats']['leaves'] for r in records for p in r['planes']))))


if __name__=='__main__':verify()
