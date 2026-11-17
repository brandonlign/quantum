import sys
if sys.flags.optimize:
    raise RuntimeError("Run without -O: this certificate uses assertions")
import re,json
from pathlib import Path
from fractions import Fraction as F
from math import comb
p=Path(__file__).resolve().parents[1]; t=(p/'manuscript/main.tex').read_text()
r=re.findall(r'\\texttt\{([IXZY]{17})\} & (\d+)',t)
assert len(r)==12
expected=[4236119552,477130815,529432071,4228805056,11524269616,10634484776,10449456558,12242979677,16154017182,11303119733,7963947440,7346381160]
enc=lambda s:sum(('IXZY'.index(a)&1)<<i|(('IXZY'.index(a)>>1)<<(17+i)) for i,a in enumerate(s))
rows=[enc(s) for s,_ in r]; assert rows==expected==[int(i) for _,i in r]
assert (p/'verification/generators.txt').read_text().splitlines()==['+'+s for s,_ in r]
m=(1<<17)-1
sp=lambda a,b:(((a&m)&(b>>17)).bit_count()+((a>>17)&(b&m)).bit_count())%2
wt=lambda a:((a&m)|(a>>17)).bit_count()
S=[0]
for a in rows:S += [b^a for b in S]
assert len(set(S))==4096 and all(sp(a,b)==0 for a in rows for b in rows)
assert wt(269748239)==5 and all(sp(269748239,a)==0 for a in rows) and 269748239 not in S
D=8065646656;s=rows[7]; assert all(sp(D,a)==0 for a in rows)
assert enc('IIIIIZYIIIIIZZZZI')==D
local=lambda a,i:((a>>i)&1)+2*((a>>(17+i))&1)
over=[i for i in range(17) if local(s,i) and local(D,i)]
assert over==[6,13,14] and all(local(s,i)==local(D,i) for i in over)
def g(a,step):
 o=0
 for i in range(17):
  x,z=(a>>i)&1,(a>>(17+i))&1
  if i<15:j=i//3*3+(i%3+step)%3
  else:j=i;x,z=z,x^z
  o|=x<<j|z<<(17+j)
 return o
print('Symmetry directions:',{step:all(g(a,step) in S for a in rows) for step in [1,2]})
y=list(map(F,['-1001/201','25/25728','101/437376','197/1749504','19/583168','5/583168','1/583168']))
C=lambda n,k:comb(n,k) if 0<=k<=n else 0
K=lambda j,i:sum((-1)**a*3**(j-a)*C(i,a)*C(17-i,j-a) for a in range(j+1))
cs=[y[0]*(i==0)+y[1]+sum(y[j+1]*(K(j,i)-4096*(i==j)) for j in range(1,6)) for i in range(18)]
nums=[0,8000,5552,3856,2337,1272,616,252,81,16,0,0,1,0,0,4,9,0]
assert cs==[F(a,3417) for a in nums] and y[0]+4096*y[1]==-1
print('PASS all rows, witness, overlap, exact bound')
# Compare the compact printed construction with its complete coordinate data.
construction=json.loads((p/'verification/normal_form/NORMAL_FORM.json').read_text())
for index,name in enumerate(['a','b']):
    line=next(line for line in t.splitlines() if line.startswith(' '+name+'&='))
    perm=list(range(17))
    cycles=re.findall(r'\(([\d\\ ]+)\)',line)
    assert cycles
    for cycle in cycles:
        cyc=list(map(int,re.findall(r'\d+',cycle)))
        for i,j in zip(cyc,cyc[1:]+cyc[:1]):perm[i]=j
    assert perm==construction['generators'][index]['permutation']
assert r"v'=\texttt{IIZXYIYXIZIIZYIII}" in t
assert enc('IIZXYIYXIZIIZYIII')==construction['seed']
assert r'14:\ X\mapsto Z\mapsto Y\mapsto X' in t
assert r'16:\ X\mapsto Y\mapsto Z\mapsto X' in t
assert all(v==[0,1,2,3] for v in construction['generators'][0]['local_label_maps'])
assert [(i,v) for i,v in enumerate(construction['generators'][1]['local_label_maps']) if v!=[0,1,2,3]]==[(14,[0,2,3,1]),(16,[0,3,1,2])]
print('PASS printed compact alternating-group construction')
