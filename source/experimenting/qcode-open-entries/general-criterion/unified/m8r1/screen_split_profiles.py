"""Exact-certificate split-enumerator screen, grouped by suffix spectrum.

The shortened C is dual to the suffix projection of S. Every suffix has
|R| stabilizer lifts, so MacWilliams gives |S| A_j(C) = sum_s K_j(wt(s_W)).
Certificates are retained as rational multipliers, not only boolean claims.
"""
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'pure'))
sys.path.insert(0,str(ROOT/'split'))
sys.path.insert(0,str(HERE.parent))
import refined_split as rs
from system import farkas_exact, verify
from model import span, split, weight


def coupled_system(n,k,d,m,minimum):
    E,b,G,V,*_=rs.build(m,minimum,n=n,k=k,d=d,pure=True)
    rows=[]
    for j in range(n-m+1):
        # The constant prefix coefficient is one for every input letter type.
        rows.append([(1 if 2*z+t==m else 2)*rs.kraw(n-m,j,w)
                     for ((z,t),w) in V])
    return np.vstack([E,np.array(rows,dtype=object)]),b,G,V


def physical_positive_control():
    # The independently reconstructed five-qubit stabilizer state retained
    # by verify_xdom_prune.py. R=<ZZZ> on the first three sites, C=0.
    S=[917,518,920,800,960]
    counts=Counter()
    for s in span(S):
        p,w=split(s,3,2)
        x=p&7;z=p>>3
        counts[((z&~x).bit_count(),x.bit_count(),weight(w,2))]+=1
    E,b,G,V=coupled_system(5,0,3,3,3)
    full_b=np.concatenate([b,np.array([32,0,0],dtype=object)])
    x=np.array([counts[(z,t,w)] for ((z,t),w) in V],dtype=object)
    assert all(v==0 for v in E@x-full_b)
    assert all(v<=0 for v in G@x)
    assert min(weight(s,5) for s in span(S) if s)==3
    return dict(status='PASS',generators=S,exact_equations=len(full_b),
                exact_inequalities=len(G),suffix_spectrum={'0':1})


def main():
    control=physical_positive_control()
    source=HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json'
    classes=json.loads(source.read_text())['classes']
    groups=defaultdict(list)
    for c in classes:
        spectrum=tuple(sorted((int(w),v) for w,v in c['representative']['weight_distribution'].items()))
        groups[spectrum].append(c['class_index'])
    E,b,G,V=coupled_system(17,3,6,8,8)
    records=[]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [source,HERE/'screen_split_profiles.py',ROOT/'pure/refined_split.py',ROOT/'split/system.py']}
    for spectrum,ids in sorted(groups.items(),key=lambda item:-len(item[1])):
        A=dict(spectrum)
        full_b=np.concatenate([b,np.array([16384*A.get(j,0) for j in range(10)],dtype=object)])
        cert,message=farkas_exact(E,full_b,G)
        record=dict(spectrum=dict(spectrum),class_indices=ids,classes=len(ids),
                    status='EXCLUDED_EXACT' if cert else 'NO_CERTIFICATE',message=str(message))
        if cert:
            nu,lam=cert
            assert verify(E,full_b,G,nu,lam)
            record['certificate']=dict(nu=list(map(str,nu)),lam=list(map(str,lam)))
        records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k not in ['class_indices','certificate']}),flush=True)
        result=dict(status='COMPLETE' if len(records)==len(groups) else 'RUNNING',
                    positive_control=control,groups_total=len(groups),records=records,
                    excluded_classes=sum(r['classes'] for r in records if r['status']=='EXCLUDED_EXACT'),
                    sha256=hashes,
                    scope='Necessary split-enumerator constraints; exact rational certificate required for exclusion')
        (HERE/'local_split_profile_screen.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
