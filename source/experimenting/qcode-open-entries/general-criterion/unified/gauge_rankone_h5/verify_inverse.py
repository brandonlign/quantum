"""Independent exhaustive inverse-label tests, including genuine distance three."""
import itertools
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
from model import embed,reduced_basis,span,weight,lift_coordinate,nullspace
from gauge_fix import fix_prefix_kernel,supported_subspace
from subcode_selection import SelectionContext
from gauge_rankone_h5.inverse import select
from mixed_frames.stage3.common import atomic_json,sha


def brute(B,n,m,h,d,cut,budget):
    ctx=SelectionContext(B,n)
    domains=[[v for v in range(1<<n) if v&((1<<h)-1)==1<<(j+1)] for j in range(h-1)]
    count=0;admissible=[]
    for rows in itertools.product(*domains):
        result=ctx.evaluate(list(rows),d,T=list(range(m)),R=[embed(((1<<m)-1)<<m,0,m,n-m)])
        if not result['distance_ok']:continue
        kept=sum(ctx.weights[c]<cut and not any((c&row).bit_count()%2 for row in rows) for c in range(1,len(ctx.words)))
        if kept<=budget:admissible.append(list(rows))
        count+=1
    return admissible,count


def main():
    report=dict(status='PASS',fixtures=[])
    saved=json.loads((ROOT/'rankone_recenter/verification.json').read_text())
    S=next(c['S'] for c in saved['positive_controls'] if c['parameters']==[7,1,3])
    fixed=fix_prefix_kernel(S,7,list(range(4)))
    assert fixed['k_fixed']==0 and fixed['local_rank_fixed']==2
    p=embed(15<<4,0,4,3);B=[p]
    for w in fixed['R_fixed']+fixed['S_fixed']:
        if len(reduced_basis(B+[w]))>len(B):B.append(w)
    fixtures=[(B,7,4,2,3)]
    for m,h,n in [(4,3,6),(5,4,7)]:
        f=n-m;R=[((1<<m)-1)<<m]
        R+=[((1<<i)|(1<<(i+1)))<<m for i in range(h-1)]
        B=[embed(r,0,m,f) for r in R]+[embed(0,2<<f,m,f),embed(15,1,m,f),embed(1<<m,1<<f,m,f)]
        assert len(B)==n and len(reduced_basis(B))==n
        fixtures.append((B,n,m,h,1))
    for B,n,m,h,d in fixtures:
        for cut,budget in [(d+1,0),(d+1,2),(d+1,1000),(n+1,(1<<(n-h+1))-2),(n+1,(1<<(n-h+1))-1)]:
            expected,count=brute(B,n,m,h,d,cut,budget)
            native=select(B,n,m,h,d,cut,budget)
            assert (native['status']=='FOUND')==bool(expected),(n,h,cut,budget,native,len(expected))
            if expected:assert native['rows'] in expected
            else:assert native['status']=='EXCLUDED_BASE'
            report['fixtures'].append(dict(n=n,m=m,h=h,d=d,cut=cut,budget=budget,admissible=len(expected),distance_admissible=count,status=native['status'],attempts=native['attempts']))
        limited=select(B,n,m,h,d,n+1,0,nodes=0)
        assert limited['status']=='LIMIT'
    report['hashes']={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'inverse.py',HERE/'inverse.cpp',HERE/'inverse',ROOT/'subcode_selection.py',ROOT/'model.py']}
    atomic_json(HERE/'inverse_verification.json',report)
    print(json.dumps(dict(status='PASS',cases=len(report['fixtures']),genuine_distance_three=True,full_label_dimensions=[1,2,3])),flush=True)


if __name__=='__main__':main()
