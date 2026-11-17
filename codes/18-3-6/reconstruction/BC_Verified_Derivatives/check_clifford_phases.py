"""Check the precise Hermitian-Pauli phases for both derivative Clifford maps."""
from pathlib import Path
import json
O=Path(__file__).parent;ns={};exec((O/'independent_binary_algebra.py').read_text(),ns);sp,solve=[ns[k] for k in ['sp','solve']]
A=json.loads((O/'BASE_CODE_A.json').read_text())['rows'];mask=(1<<18)-1
if not __debug__:raise RuntimeError('Assertions required')
def phase_product(a,b,n):
 m=(1<<n)-1;qa=((a&m)&(a>>n)).bit_count();qb=((b&m)&(b>>n)).bit_count();c=a^b;qc=((c&m)&(c>>n)).bit_count()
 return (qa+qb-qc+2*((a>>n)&(b&m)).bit_count())%4
# Independent explicit one-qubit matrices verify all sixteen Pauli products.
I=((1,0),(0,1));X=((0,1),(1,0));Z=((1,0),(0,-1));Y=((0,-1j),(1j,0));P=[I,X,Z,Y]
def mul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)) for i in range(2))
for a in range(4):
 for b in range(4):
  phase=1j**phase_product(a,b,1);assert mul(P[a],P[b])==tuple(tuple(phase*x for x in row) for row in P[a^b])
records=[]
for name in ['code_b_style','code_c_style']:
 D=O/name;T=json.loads((D/'TRANSFORM.json').read_text());v=T['v'];R=T['rows'];signs=[]
 for r,t in zip(A,R):
  if not sp(v,r,18):assert t==r;signs.append(0)
  else:
   assert t==r^v;phase=(1+phase_product(v,r,18))%4;assert phase in [0,2];signs.append(phase//2)
 eq=[(r>>18)|((r&mask)<<18) for r in R];w=solve(eq,signs,36)[0]
 assert all(sp(w,r,18)==s for r,s in zip(R,signs))
 records.append({'code':name,'transvection_v':v,'unitary_convention':'U_v = exp(+i*pi*P_v/4) = (I+i P_v)/sqrt(2)','bare_conjugation_negative_generator_indices_zero_based':[i for i,s in enumerate(signs) if s],'Pauli_correction_w':w,'Pauli_correction_string':''.join('IXZY'[(w>>i&1)+2*(w>>(18+i)&1)] for i in range(18)),'exact_relation':'(P_w U_v) P_(A_j) (P_w U_v)^dagger = P_(derivative_j), with every displayed generator Hermitian and positive'})
(O/'PHASE_CORRECTIONS.json').write_text(json.dumps({'status':'PASS','single_qubit_Pauli_product_matrix_checks':16,'global_generator_phase_checks':30,'derivatives':records},indent=2)+'\n');print(json.dumps(records,indent=2))
