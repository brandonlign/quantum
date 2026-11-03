"""Physical positive control for audit16's independently generated input."""
from audit16 import OUT, PKG, nullspace, quotient_distances, physical_constraints, search_binary
import json
import subprocess

assert __debug__, 'Run without -O'
OUT.mkdir(exist_ok=True)
build, binary = search_binary()

n=16
mask=(1<<n)-1
def weight(w):
    return ((w & mask)|(w>>n)).bit_count()
def span(gens):
    words=[0]
    for g in gens:
        words += [w^g for w in words]
    return words
gens=[]
for line in (PKG/'search/q1645.txt').read_text().splitlines():
    x,z=line.split('|')
    gens.append(sum(int(b)<<i for i,b in enumerate(x.split())) | sum(int(b)<<(n+i) for i,b in enumerate(z.split())))
assert nullspace(gens,32)[0]==12
equations=[(w>>n)|((w&mask)<<n) for w in gens]
assert all((a&b).bit_count()%2==0 for a in equations for b in gens)
S=set(span(gens))
dual_gens=nullspace(equations,32)[1]
assert len(dual_gens)==20
logical_distance=min(weight(w) for w in span(dual_gens) if w not in S)
assert logical_distance==5
p=min(S-{0},key=lambda w:(weight(w),w))
m=weight(p)
assert m==6
support=[i for i in range(n) if ((p&mask)|(p>>n))>>i&1]
rest=[i for i in range(n) if i not in support]
f=len(rest)
projection=[]
for w in gens:
    x=sum(((w>>i)&1)<<j for j,i in enumerate(rest))
    z=sum(((w>>(i+n))&1)<<j for j,i in enumerate(rest))
    projection.append(z|(x<<f))
C=nullspace(projection,2*f)[1]
assert len(C)==9
dist=quotient_distances(C,f)
assert len(dist)==2048
rows=physical_constraints(m)
data=f'{m} {len(dist)}\n'+' '.join(map(str,dist))+f'\n{len(rows)}\n'+''.join(f'{c} {t}\n' for c,t in rows)
result=subprocess.run([str(binary)],input=data,text=True,capture_output=True,timeout=120)
assert result.returncode==0 and result.stdout.startswith('FOUND'),result.stdout
report={'status':'PASS','physical_parameters':[16,4,5],'stabilizer_rank':12,'centralizer_rank':20,'minimum_stabilizer_weight':m,'suffix_dimension':len(C),'search':result.stdout.strip(),'scope':'Independent Pauli-space reconstruction and quotient BFS; no project proof modules imported.'}
(OUT/'control16_independent.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
build.cleanup()
