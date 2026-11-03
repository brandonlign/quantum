"""Complete one-coordinate lifts from the existing length8,dim6,d5 census.

Outputs are extension representatives, NOT inequivalent length9 classes.
"""
from collections import Counter
import hashlib,json,time
from pathlib import Path
import numpy as np
from model import gram,reduced_basis,weight,span

BASE=Path(__file__).resolve().parent
DATA=BASE.parents[1]/'release16/census/cells.json'


def maps(c):
    out=[(0,0)]+[(a,0) for a in range(1,1<<c)]
    out += [(a,b) for a in range(1,1<<c) for b in range(a+1,1<<c) if b<a^b]
    return out


def lift(gens,a,b,f):
    mask=(1<<f)-1
    return [(g&mask)|(((a>>i)&1)<<f)|
            (((g>>f)|(((b>>i)&1)<<f))<<(f+1)) for i,g in enumerate(gens)]


def main():
    started=time.monotonic();parents=json.loads(DATA.read_text())['8,6'];forms=maps(6)
    assert len(forms)==715 and len(set(forms))==715
    parity=np.array([[(a&w).bit_count()%2 for w in range(64)] for a in range(64)],dtype=bool)
    coverage=np.asarray([parity[a]|parity[b] for a,b in forms])
    hist=Counter();counts=[];lifts=[]
    for i,gens in enumerate(parents):
        assert len(gens)==len(reduced_basis(gens))==6
        ws=span(gens);assert len(ws)==64 and min(weight(w,8) for w in ws if w)>=5
        # Recompute indexed words explicitly: coefficient labels match map rows.
        indexed=[0]
        for g in gens:indexed += [w^g for w in indexed]
        bad=[j for j,w in enumerate(indexed) if weight(w,8)==5]
        valid=np.flatnonzero(coverage[:,bad].all(axis=1))
        local=Counter()
        for idx in valid:
            a,b=forms[int(idx)];new=lift(gens,a,b,8)
            assert len(reduced_basis(new))==6
            actual=min(weight(w,9) for w in span(new) if w)
            assert actual>=6
            h=6-len(reduced_basis(gram(new,9)))
            local[h]+=1;hist[h]+=1
            lifts.append(dict(parent=i,tail_forms=[a,b],generators=new,hull=h,distance=actual))
        counts.append(dict(parent=i,weight5_words=len(bad),extensions=len(valid),hull_counts=dict(local)))
        if i%200==0:print(json.dumps(dict(status='RUNNING',parents=i+1,extensions=len(lifts),hull_counts=dict(hist))),flush=True)
    output=dict(status='COMPLETE',parent_cell='8,6',parent_count=len(parents),
                canonical_maps_per_parent=len(forms),extension_representatives=len(lifts),
                hull_counts=dict(hist),counts=counts,lifts=lifts,
                target_branch='[[17,3,6]],m=8,r=2',required_hulls=[0,2,4],
                branch_excluded=not any(hist[h] for h in [0,2,4]),
                seconds=time.monotonic()-started,
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [DATA,Path(__file__),BASE/'model.py']},
                scope='Complete lifts relative to the supplied complete parent census; '
                      'representatives may be equivalent. Does not settle other local branches.')
    (BASE/'suffix9_extensions.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['counts','lifts','sha256']}))
if __name__=='__main__':main()
