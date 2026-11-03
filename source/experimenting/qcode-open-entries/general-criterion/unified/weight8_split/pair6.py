"""Exact Farkas: pure [[17,3,6]] with A6=2, A7=0 (only allowed A6=2 profile). Two weight-6 words
w1=Z^6, w2 overlapping t sites, s of them Z (cancel), t-s X; wt(w1+w2)=12-t-s>=8. Union block."""
import sys, json, time; sys.path.insert(0,'split')
from system import build, farkas_exact
n,k,d=17,3,6
out={}
for t,s in [(4,0),(3,1),(2,2),(2,0),(1,1),(0,0)]:
    u=12-t
    w1='Z'*6+'I'*(u-6)
    # sites 0..5 = w1; overlap = sites 6-t..5 ; w2-only = 6..u-1
    w2=['I']*u
    for i in range(6-t,6-t+s): w2[i]='Z'
    for i in range(6-t+s,6): w2[i]='X'
    for i in range(6,u): w2[i]='Z'
    w2=''.join(w2); t0=time.time()
    Meq,b,G,cols,*_=build(n,k,d,[(u,[w1,w2])],weights={1:0,2:0,3:0,4:0,5:0,6:2,7:0})
    cert,msg=farkas_exact(Meq,b,G)
    out[f't{t}s{s}']=dict(w1=w1,w2=w2,cert=bool(cert),cols=len(cols),sec=round(time.time()-t0))
    print(out[f't{t}s{s}'],flush=True)
    json.dump(out,open('general-criterion/unified/weight8_split/pair6.json','w'),indent=1)
