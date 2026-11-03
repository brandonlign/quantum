"""Necessary capacity screen for branch m8,r1 after gauge-fixing ell=(hull+1)/2 prefix logicals.

Fixed problem has zero kernel, so M: E+ -> Q=im Gamma is a bijection and
#{0!=s in Q: delta(s)<=j} <= #{0!=e in E+: cost+(e)>=6-j}, j=0..5.
A class fails for a profile if any j violates this; it is EXCLUDED by this screen only if it
fails for EVERY admissible fixed-logical profile. Profiles: fixed_profiles_12.json.
"""
import json,sys
import numpy as np
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from model import gram,span
T=np.load(HERE/'delta_tables.npy');meta=json.loads((HERE/'delta_meta.json').read_text())
prof=json.loads((HERE/'fixed_profiles_12.json').read_text()) if (HERE/'fixed_profiles_12.json').exists() else {}
S=json.loads(Path(meta['source']).read_text())['classes']
ell=int(sys.argv[1]);profiles=[tuple(p) for p in prof.get(str(ell),[])]
if not profiles:sys.exit('no profiles for ell=%d'%ell)
res=dict(ell=ell,profiles=profiles,excluded=[],survivors=[],slack={})
for ci,c in enumerate(S):
    h=c['representative']['hull']
    if (h+1)//2!=ell:continue
    G=gram(meta['P_bases'][ci],9);Q=np.array(span(G));assert len(Q)==1<<(13-h)
    dq=T[ci][Q[Q!=0]]
    act=[int((dq<=j).sum()) for j in range(6)]
    ok=[p for p in profiles if all(act[j]<=p[6-j] for j in range(6))]
    (res['survivors'] if ok else res['excluded']).append(ci)
    res['slack'][ci]=act
print(json.dumps(dict(ell=ell,classes=len(res['slack']),excluded=len(res['excluded']),survivors=len(res['survivors']))))
ex=next(iter(res['slack'].values()));print('sample actual counts',ex,'profiles',profiles)
json.dump(res,open(HERE/f'capacity_ell{ell}.json','w'))
