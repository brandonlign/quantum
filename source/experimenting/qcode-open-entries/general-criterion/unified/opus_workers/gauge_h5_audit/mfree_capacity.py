"""Bounded audit fixture: M-free heavy-coset capacity test for (H,p,C).

For a heavy S_plus coset with H-signature sigma (bit0=0, sigma!=0) and
C-signature tau, every bijection e -> s_tau + M e must satisfy
cost_sigma(e) + delta(s) >= 6. Counting gives, for j=0..5,
  #{H-cosets in u_sigma+H^perp of weight <= j} <= #{s in s_tau+im Gamma: delta(s) >= 6-j}.
Necessary only. Reads the three pilot (R,C) pairs; no search.
"""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from model import ExtensionModel,span,symp,weight,reduced_basis,lift_coordinate
from mixed_frames.frames import adapted_frame

def prefix_classes(R,m=8):
    p=255<<8;B=[p]
    for r in reduced_basis(R):
        if len(reduced_basis(B+[r]))>len(B):B.append(r)
    Hs=span(B);out={}
    for u in range(1<<(2*m)):
        sig=sum(symp(b,u,m)<<i for i,b in enumerate(B))
        if sig&1 or sig==0:continue
        key=min(u^h for h in Hs)
        if key==u:
            out.setdefault(sig,[]).append(min(weight(u^h,m) for h in Hs))
    return out

def main():
    rows=[]
    for i in range(3):
        d=json.loads((ROOT/f'gauge_rankone_h5/pilot_{i}.json').read_text())
        R,C=d['R'],d['C'];frame=adapted_frame(8,R)
        model=ExtensionModel(17,0,6,8,frame.R,C)
        im=set(span(model.Gamma));assert len(im)==256
        reps=[];seen=set()
        for s in range(1<<13):
            if s not in seen:reps.append(s);seen|={s^t for t in im}
        assert len(reps)==32
        supply=[[sum(model.distances[r^t]>=6-j for t in im) for j in range(6)] for r in reps]
        cls=prefix_classes(frame.R)
        res={}
        for sig,ws in sorted(cls.items()):
            assert len(ws)==256
            demand=[sum(w<=j for w in ws) for j in range(6)]
            ok=[t for t,sup in enumerate(supply) if all(a<=b for a,b in zip(demand,sup))]
            res[sig]=dict(demand=demand,passing_tau=len(ok))
        rows.append(dict(pilot=i,signatures=res,
                         excluded_by_mfree=any(v['passing_tau']==0 for v in res.values())))
    print(json.dumps(rows,indent=1))
    (Path(__file__).with_name('mfree_capacity.json')).write_text(json.dumps(rows,indent=1)+'\n')
main()
