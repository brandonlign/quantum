from pathlib import Path
import json,sys
O=Path(__file__).parent;ns={};exec((O/'independent_binary_algebra.py').read_text(),ns)
sp,wt,basis,perp,span=[ns[k] for k in ['sp','wt','basis','perp','span']]
A=json.loads((O/'BASE_CODE_A.json').read_text())['rows'];mask=int(sys.argv[1]) if len(sys.argv)>1 else sum(1<<(18+i) for i in [2,3,4,5]);name=sys.argv[2] if len(sys.argv)>2 else 'specified_z2345';D=O/name;D.mkdir(exist_ok=True)
S=[r^(mask if sp(r,mask,18) else 0) for r in A];assert len(basis(S))==15 and all(not sp(a,b,18) for a in S for b in S);N=perp(S,18);assert len(N)==21
(D/'ROWS.txt').write_text('\n'.join(map(str,S))+'\n');(D/'INPUT.txt').write_text(' '.join(map(str,S+N))+'\n');(D/'TRANSFORM.json').write_text(json.dumps({'base':'Code A','operation':'binary symplectic transvection r -> r + sp(r,v)*v','v':mask,'Pauli_v':''.join('IXZY'[(mask>>i&1)+2*(mask>>(18+i)&1)] for i in range(18)),'rows':S},indent=2)+'\n')
print(D)
