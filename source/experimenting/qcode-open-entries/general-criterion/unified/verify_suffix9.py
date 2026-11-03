"""Unrestricted-map verification independent of census canonicalization/algebra."""
from pathlib import Path
import json,hashlib,time
from collections import Counter
import numpy as np
BASE=Path(__file__).resolve().parent
DATA=BASE.parents[1]/'release16/census/cells.json'

def rank(rows):
    pivots={}
    for row in rows:
        while row:
            bit=row.bit_length()-1
            if bit in pivots:row^=pivots[bit]
            else:pivots[bit]=row;break
    return len(pivots)

def sp(a,b,n):
    mask=(1<<n)-1
    return (((a&mask)&(b>>n)).bit_count()+((b&mask)&(a>>n)).bit_count())%2

def physical_weight(w,n):return ((w&((1<<n)-1))|(w>>n)).bit_count()

def canonical(a,b):
    vals=sorted({a,b,a^b}-{0})
    return (0,0) if not vals else (vals[0],0) if len(vals)==1 else tuple(vals[:2])

def main():
    t=time.monotonic();parents=json.loads(DATA.read_text())['8,6'];saved=json.loads((BASE/'suffix9_extensions.json').read_text())
    actual={}
    for item in saved['lifts']:actual.setdefault(item['parent'],{})[tuple(item['tail_forms'])]=item
    forms=[(a,b) for a in range(64) for b in range(64)]
    parity=np.array([[(a&w).bit_count()%2 for w in range(64)] for a in range(64)],dtype=bool)
    coverage=np.asarray([parity[a]|parity[b] for a,b in forms])
    unrestricted=0;hulls=Counter();heavy_survivors=[];checked=0
    for idx,gens in enumerate(parents):
        words=[0]
        for g in gens:words += [w^g for w in words]
        assert len(set(words))==64
        bad=[j for j,w in enumerate(words) if physical_weight(w,8)==5]
        good=np.flatnonzero(coverage[:,bad].all(axis=1));unrestricted+=len(good)
        orbit=Counter(canonical(*forms[int(i)]) for i in good)
        assert set(orbit)==set(actual.get(idx,{})),idx
        for tail,multiplicity in orbit.items():
            a,b=tail;assert multiplicity==(1 if a==b==0 else 3 if b==0 else 6)
            new=[(g&255)|(((a>>i)&1)<<8)|((g>>8)<<9)|(((b>>i)&1)<<17) for i,g in enumerate(gens)]
            assert new==actual[idx][tail]['generators']
            expanded=[0]
            for g in new:expanded += [w^g for w in expanded]
            assert len(set(expanded))==64 and min(physical_weight(w,9) for w in expanded if w)>=6
            G=[sum(sp(u,v,9)<<j for j,v in enumerate(new)) for u in new]
            h=6-rank(G);assert h==actual[idx][tail]['hull'];hulls[h]+=1
            rad=[w for w in expanded if all(sp(w,g,9)==0 for g in new)]
            assert len(rad)==1<<h
            if h==2:assert sorted(physical_weight(w,9) for w in rad if w)==[6,6,6]
            if h==0:heavy_survivors.append(actual[idx][tail])
            checked+=1
    assert checked==284 and hulls=={0:98,2:171,6:15}
    # The target has q=h/2 and dim N=q; N<=rad C, minimum weight8.
    # h6 contradicts Plotkin's |N|<=floor(32/5)=6; h4 absent;
    # h2 has no permissible nonzero N word.
    report=dict(status='PASS',parent_classes=len(parents),all_tail_maps_per_parent=4096,
        total_tail_maps_checked=len(parents)*4096,valid_unrestricted_maps=unrestricted,
        canonical_representatives=checked,hull_counts=dict(hulls),
        remaining_target_representatives=len(heavy_survivors),remaining_hull=0,
        remaining_suffix_stabilizer_dimension=0,
        scope='Complete relative to the supplied complete parent census; '
              'independent unrestricted tail-map and physical radical checks, '
              'not a full stabilizer extension or new table-entry exclusion.',
        seconds=time.monotonic()-t,
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [Path(__file__),DATA,BASE/'suffix9_extensions.json']})
    (BASE/'suffix9_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    (BASE/'suffix9_target_candidates.json').write_text(json.dumps(dict(status='COMPLETE_REDUCTION',
        target='[[17,3,6]],m=8,r=2',representatives=heavy_survivors),indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
