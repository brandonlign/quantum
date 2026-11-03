"""Letter-refined split LP for pure [[17,3,6]] with Z^8 in S, min stabilizer weight 8,
plus the exact weight distribution of the shortened code C = {suffix-only words of S^perp}:
B(0,0,b) = A_b(C). One LP per distinct (weight distribution) among the 32,865 classes.
CERT = exact Farkas certificate verified in rationals (refined_split.farkas_exact)."""
import sys,json
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'pure'));sys.path.insert(0,str(ROOT/'split'))
import refined_split as rs
from system import farkas_exact
n,k,d,m=17,3,6,8;f=n-m;MS=2**(n-k)
def build_C(Aw):
    Meq,b,G,V,eqlab,lelab=rs.build(m,m,n=n,k=k,d=d,pure=True)
    BK={(z,t):rs.pmul(rs.pmul(rs.ppow(rs.PI,m-z-t),rs.ppow(rs.PZ,z)),rs.ppow(rs.PX,t)) for z in range(m+1) for t in range(m+1-z)}
    rows=[];rhs=[]
    for bc in range(f+1):
        r=[0]*len(V)
        for j,((key,bb)) in enumerate(V):
            zs={(key[0],key[1]),(m-key[0]-key[1],key[1])}
            for (z,t) in zs:r[j]+=BK[(z,t)].get((0,0),0)*rs.kraw(f,bc,bb)
        rows.append(r);rhs.append(MS*Aw.get(bc,0))
    Meq=np.vstack([Meq,np.array(rows,dtype=object)]);b=np.concatenate([b,np.array(rhs,dtype=object)])
    return Meq,b,G
if __name__=='__main__':
    S=json.load(open(ROOT/'general-criterion/unified/opus_parallel/suffix9_dimension5/suffix9_dim5.json'))['classes']
    groups=defaultdict(list)
    for c in S:groups[tuple(sorted((int(a),v) for a,v in c['representative']['weight_distribution'].items()))].append(c['class_index'])
    # control: no C constraint
    Meq,b,G,*_=rs.build(m,m,n=n,k=k,d=d,pure=True);c0,_=farkas_exact(Meq,b,G);print('baseline m8 (no C):','CERT' if c0 else 'none',flush=True)
    out=[]
    for wd,ids in sorted(groups.items(),key=lambda x:-len(x[1])):
        Meq,b,G=build_C(dict(wd));cert,msg=farkas_exact(Meq,b,G)
        st=rs.lp_status(Meq,b,G)
        out.append(dict(weights=wd,classes=len(ids),cert=bool(cert),lp_status=int(st)))
        print(len(ids),wd,'CERT' if cert else 'none','lp',st,flush=True)
    json.dump(dict(results=out,certified_classes=sum(o['classes'] for o in out if o['cert'])),open('split_with_C.json','w'),indent=1)
    print('certified classes',sum(o['classes'] for o in out if o['cert']),'of',len(S))
