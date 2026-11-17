"""Portable check of the local physical tables and reconstruction of CodeA."""
from pathlib import Path
import json
p=Path(__file__).parent;ctx=json.loads((p/'PHYSICAL_HIT_CONTEXT.json').read_text());F=json.loads((p/'HIT.json').read_text())['map'];L,R=ctx['left'],ctx['right']
def span(rows):
 w=[0]
 for v in rows:w +=[x^v for x in w]
 return w
def wt(v,n):return ((v&((1<<n)-1))|(v>>n)).bit_count()
def embed(v,n,off):return ((v&((1<<n)-1))<<off)|((v>>n)<<(18+off))
def kernel(rows,width):
 piv={}
 for v in rows:
  while v:
   i=v.bit_length()-1
   if i in piv:v^=piv[i]
   else:piv[i]=v;break
 out=[]
 for j in range(width):
  if j in piv:continue
  v=1<<j
  for i in sorted(piv):
   if (piv[i]&v).bit_count()%2:v^=1<<i
  out.append(v)
 return out
for local in [L,R]:
 n=local['n'];cw=span(local['C']);qw=span(local['Q']);assert len(set(cw))==len(cw);assert all(wt(v,n)>=6 for v in cw[1:])
 assert [min(wt(q^c,n) for c in cw) for q in qw]==local['cost']
assert sorted(F)==list(range(1024));assert F==span([F[1<<i] for i in range(10)])
assert all(L['cost'][v]+R['cost'][F[v]]>=6 for v in range(1,1024))
rq=span(R['Q']);N=[embed(v,8,0) for v in L['C']]+[embed(v,10,8) for v in R['C']]+[embed(L['Q'][i],8,0)^embed(rq[F[1<<i]],10,8) for i in range(10)]
mask=(1<<18)-1;S=kernel([(v>>18)|((v&mask)<<18) for v in N],36);expected=list(map(int,(p.parent/'RAW15.rows').read_text().split()));assert expected[:2]==[18,15] and S==expected[2:]
print(json.dumps({'local_coset_tables_checked':True,'all_1023_nonzero_cost_inequalities_pass':True,'map_linear_and_bijective':True,'reconstructed_stabilizer_rows_match':True,'raw15':S},indent=2))
