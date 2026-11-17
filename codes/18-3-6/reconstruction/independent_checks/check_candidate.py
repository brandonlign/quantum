from pathlib import Path
import json,hashlib,collections
O=Path(__file__).parent;ns={};exec((O/'independent_binary_algebra.py').read_text(),ns)
basis,span,sp,wt,perp=[ns[k] for k in ['basis','span','sp','wt','perp']]
S=[26277347457,59444428930,52902924420,3136815240,10548711568,27939967136,16496890048,63874044160,21506327040,66477130752,15987673088,20945838080,51785547776,16911630336,66846720]
assert len(S)==15 and all(0<x<1<<36 for x in S) and len(basis(S))==15
assert all(sp(a,b,18)==0 for a in S for b in S)
N=perp(S,18);assert len(N)==len(basis(N))==21
assert all(not sp(a,b,18) for a in S for b in N)
W=set(span(S));assert len(W)==32768
A=collections.Counter(wt(x,18) for x in W)
(O/'RAW15.rows').write_text('\n'.join(map(str,S))+'\n')
(O/'NORMALIZER.rows').write_text('\n'.join(map(str,N))+'\n')
(O/'INPUT.txt').write_text(' '.join(map(str,S+N))+'\n')
J={'rank':15,'isotropic':True,'normalizer_rank':21,'stabilizer_minimum':min(w for w in A if w),'stabilizer_weight_enumerator':dict(sorted(A.items())),'raw_rows_sha256':hashlib.sha256((O/'RAW15.rows').read_bytes()).hexdigest()}
(O/'PRELIMINARY.json').write_text(json.dumps(J,indent=2)+'\n');print(json.dumps(J,indent=2))
