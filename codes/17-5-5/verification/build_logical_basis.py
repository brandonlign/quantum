"""Construct and verify a symplectic logical Pauli basis; no external packages."""
from pathlib import Path
import json
P=Path(__file__).parent
rows=json.loads((P/'generators.json').read_text())['rows']
n=17;mask=(1<<n)-1
def bracket(a,b):
    return (((a&mask)&(b>>n)).bit_count()+((a>>n)&(b&mask)).bit_count())%2
def insert(v,basis):
    while v:
        i=v.bit_length()-1
        if i in basis:v^=basis[i]
        else:basis[i]=v;return True
    return False
constraints={}
for s in rows:insert((s>>n)|((s&mask)<<n),constraints)
normalizer=[]
for j in range(2*n):
    if j in constraints:continue
    v=1<<j
    for i in sorted(constraints):
        if (constraints[i]&v).bit_count()%2:v^=1<<i
    normalizer.append(v)
span={}
for s in rows:insert(s,span)
logical=[]
for v in normalizer:
    if insert(v,span):logical.append(v)
if len(logical)!=10:raise RuntimeError('quotient dimension')
pairs=[]
while logical:
    a=logical.pop(0)
    j=next(i for i,b in enumerate(logical) if bracket(a,b))
    b=logical.pop(j)
    logical=[v^(a if bracket(v,b) else 0)^(b if bracket(v,a) else 0) for v in logical]
    pairs.append((a,b))
flat=[v for pair in pairs for v in pair]
if not all(bracket(v,s)==0 for v in flat for s in rows):raise RuntimeError('centralizer')
if not all(bracket(a,b)==int(i//2==j//2 and i!=j) for i,a in enumerate(flat) for j,b in enumerate(flat)):raise RuntimeError('canonical brackets')
span={}
for v in rows+flat:insert(v,span)
if len(span)!=22:raise RuntimeError('independence modulo stabilizer')
def describe(v):
    return {'integer':v,'pauli':''.join('IXZY'[((v>>i)&1)+2*((v>>(17+i))&1)] for i in range(n)), 'weight':((v&mask)|(v>>n)).bit_count()}
out={'status':'PASS','convention':'site 0 first; X bits 0..16, Z bits 17..33', 'logical_pairs':[{'X':describe(a),'Z':describe(b)} for a,b in pairs], 'scope':'One valid logical basis, not weight-minimized or a claimed optimal encoder.'}
(P/'logical_basis.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS: five canonical logical X/Z pairs, rank 22 with stabilizer.')
