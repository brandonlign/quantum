"""Independent unrestricted tail-map replay of Opus's length11 lift census."""
from pathlib import Path
from collections import Counter
import hashlib,json,time
import numpy as np
BASE=Path(__file__).resolve().parent
DATA=BASE.parents[1]/'release16/census/cells.json'
OPUS=BASE/'opus_parallel/parent_extension_11/enumeration.json'

def rank(rows):
    piv={}
    for v in rows:
        while v:
            bit=v.bit_length()-1
            if bit in piv:v^=piv[bit]
            else:piv[bit]=v;break
    return len(piv)

def wt(v,n):return ((v&((1<<n)-1))|(v>>n)).bit_count()

def sp(u,v,n):
    mask=(1<<n)-1
    return (((u&mask)&(v>>n)).bit_count()+((v&mask)&(u>>n)).bit_count())%2

def canonical(a,b):
    vals=sorted({a,b,a^b}-{0})
    return (0,0) if not vals else (vals[0],0) if len(vals)==1 else tuple(vals[:2])

def main():
    t=time.monotonic();source=json.loads(OPUS.read_text());parents=json.loads(DATA.read_text())['10,10']
    assert len(source['per_parent'])==len(parents)==37
    assert source['source_sha256']==hashlib.sha256(DATA.read_bytes()).hexdigest()
    expected=set(tuple(v) for v in source['survivors']);found=set();unrestricted=0;records=[]
    for pi,g in enumerate(parents):
        words=[0]
        for v in g:words += [w^v for w in words]
        assert len(set(words))==1024 and min(wt(w,10) for w in words if w)==5
        bad=[i for i,w in enumerate(words) if wt(w,10)==5]
        assert len(bad)<=128
        zero_lo=np.zeros(1024,dtype=np.uint64);zero_hi=zero_lo.copy()
        for a in range(1024):
            lo=hi=0
            for j,w in enumerate(bad):
                if (a&w).bit_count()%2==0:
                    if j<64:lo |= 1<<j
                    else:hi |= 1<<(j-64)
            zero_lo[a]=lo;zero_hi[a]=hi
        valid=((zero_lo[:,None]&zero_lo[None,:])==0)&((zero_hi[:,None]&zero_hi[None,:])==0)
        ii,jj=np.nonzero(valid);unrestricted+=len(ii)
        orbits=Counter(canonical(int(a),int(b)) for a,b in zip(ii,jj))
        for (a,b),count in orbits.items():
            assert count==(1 if a==b==0 else 3 if b==0 else 6)
            found.add((pi,a,b));new=[(w&1023)|(((a>>i)&1)<<10)|((w>>10)<<11)|(((b>>i)&1)<<21) for i,w in enumerate(g)]
            expanded=[0]
            for v in new:expanded += [w^v for w in expanded]
            assert len(set(expanded))==1024 and min(wt(w,11) for w in expanded if w)==6
            G=[sum(sp(u,v,11)<<j for j,v in enumerate(new)) for u in new]
            h=10-rank(G);records.append(dict(parent=pi,tail_forms=[a,b],generators=new,
                hull=h,distance=6,even=all(wt(w,11)%2==0 for w in expanded)))
    assert found==expected and len(found)==3
    assert sorted(r['hull'] for r in records)==[0,10,10]
    report=dict(status='PASS',parent_classes=37,unrestricted_maps_per_parent=1048576,
        total_unrestricted_maps_checked=37*1048576,valid_unrestricted_maps=unrestricted,
        canonical_extension_representatives=len(records),required_hull=4,
        hull_counts=dict(Counter(r['hull'] for r in records)),
        target_branch='[[17,3,6]],m=6,r=2',branch_excluded=True,records=records,
        seconds=time.monotonic()-t,
        scope='Independent unrestricted tail-map enumeration and physical weight/Gram checks. '
              'Complete relative to the supplied complete37-parent census; excludes one '
              'whole local branch, not the complete parameter set.',
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),DATA,OPUS]})
    (BASE/'suffix11_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['records','sha256']}))
if __name__=='__main__':main()
