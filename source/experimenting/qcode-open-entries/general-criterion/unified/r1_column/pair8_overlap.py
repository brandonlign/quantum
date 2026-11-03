"""Exact Farkas screen, rank-one weight-eight column of pure [[17,3,6]].
A8>=120 words of weight 8 on 17 sites: site degrees sum to >=960, so the mean pairwise support
overlap is >=26628/7140>3; some pair p,q has overlap t>=4. t=8 would give local rank 2 (excluded
in this column). Normalize p=Z^8; q has s Z-letters and t-s X-letters on the overlap (t-s even),
Z elsewhere; wt(p+q)=16-t-s>=6. Block = union (16-t sites) with R=<p,q>, fixed profile (A6,A7).
CERT = the necessary relaxation is infeasible for that geometry/profile."""
import sys,json,time;sys.path.insert(0,'split')
from system import build,farkas_exact
n,k,d=17,3,6
geoms=[(t,s) for t in range(4,8) for s in range(t%2,t+1,2) if 16-t-s>=6]
profiles=[(0,0),(0,1),(0,2),(1,0),(1,1),(2,0)]
out={}
for t,s in geoms:
    u=16-t
    p='Z'*8+'I'*(8-t)
    q=['I']*(8-t)+['Z']*s+['X']*(t-s)+['Z']*(8-t)
    q=''.join(q)
    assert len(p)==len(q)==u
    for a6,a7 in profiles:
        w=16-t-s
        if w==6 and a6==0 or w==7 and a7==0: out[f't{t}s{s}_{a6}{a7}']='IMPOSSIBLE_BY_COUNT';continue
        t0=time.time()
        A,b,G,cols,*_=build(n,k,d,[(u,[p,q])],weights={1:0,2:0,3:0,4:0,5:0,6:a6,7:a7})
        cert,msg=farkas_exact(A,b,G)
        out[f't{t}s{s}_{a6}{a7}']='CERT' if cert else 'open'
        print(t,s,a6,a7,out[f't{t}s{s}_{a6}{a7}'],len(cols),round(time.time()-t0,1),flush=True)
        json.dump(out,open('general-criterion/unified/r1_column/pair8_overlap.json','w'),indent=1)
