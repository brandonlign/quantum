"""Self-contained binary stabilizer construction; Python 3.10+, standard library.
Coordinates in the displayed cycles are numbered 1..17. Pauli phases omitted.
"""
from itertools import combinations, product
from collections import Counter
import json
if not __debug__:
 raise RuntimeError("Run this verifier without Python optimization flags.")
N=17; MASK=(1<<N)-1

def permutation(cycles):
 p=list(range(N))
 for cyc in cycles:
  for i,j in zip(cyc,cyc[1:]+cyc[:1]):p[i-1]=j-1
 return p
A=permutation([[1,11],[2,10],[3,12],[4,7],[5,9],[6,8],[14,15],[16,17]])
B=permutation([[1,6,10],[2,12,8],[3,7,4],[5,9,11],[13,14,16]])
def act(v,p,twists={}):
 out=0
 for i in range(N):
  t=((v>>i)&1)|(((v>>(N+i))&1)<<1)
  if i in twists:t=twists[i][t]
  out|=(t&1)<<p[i]|(t>>1)<<(N+p[i])
 return out
def a(v):return act(v,A)
def b(v):return act(v,B,{14:(0,2,3,1),16:(0,3,1,2)})
def parse(word):
 assert len(word)==N
 return sum(('IXZY'.index(t)&1)<<i|('IXZY'.index(t)>>1)<<(N+i) for i,t in enumerate(word))
def weight(v):return ((v&MASK)|(v>>N)).bit_count()
def sp(v,w):return (((v&MASK)&(w>>N)).bit_count()+((w&MASK)&(v>>N)).bit_count())%2
seed=parse('IIZXYIYXIZIIZYIII')
orb={seed};queue=[seed]
for v in queue:
 for op in [a,b]:
  w=op(v)
  if w not in orb:orb.add(w);queue.append(w)
assert len(orb)==30 and all(weight(v)==8 for v in orb)
assert all(sp(v,w)==0 for v in orb for w in orb)
rows=[];echelon={}
for x in sorted(orb):
 y=x
 while y:
  k=y.bit_length()-1
  if k in echelon:y^=echelon[k]
  else:echelon[k]=y;rows.append(x);break
assert len(rows)==12
S={0}
for r in rows:S|={v^r for v in S}
assert len(S)==4096
swe=dict(sorted(Counter(weight(v) for v in S).items()))
assert swe=={0:1,8:30,10:720,12:1344,14:1680,16:321}
def syndrome(v):return sum(sp(v,r)<<i for i,r in enumerate(rows))
seen={syndrome(0)}
for wt in [1,2]:
 for sites in combinations(range(N),wt):
  for labels in product(range(1,4),repeat=wt):
   v=sum((t&1)<<i|(t>>1)<<(N+i) for i,t in zip(sites,labels))
   h=syndrome(v);assert h not in seen;seen.add(h)
assert len(seen)==1276  # Injectivity proves no nonzero normalizer word of wt<=4.
witness=None
for sites in combinations(range(N),5):
 for labels in product(range(1,4),repeat=5):
  v=sum((t&1)<<i|(t>>1)<<(N+i) for i,t in zip(sites,labels))
  if syndrome(v)==0:
   assert v not in S;witness=v;break
 if witness is not None:break
assert witness is not None
# Group relations checked on a basis of the full ambient binary space.
for v in [1<<i for i in range(2*N)]:
 assert a(a(v))==v and b(b(b(v)))==v
 w=v
 for _ in range(5):w=a(b(w))
 assert w==v
identity=tuple(1<<i for i in range(2*N));group={identity};todo=[identity]
for h in todo:
 for op in [a,b]:
  g=tuple(op(v) for v in h)
  if g not in group:group.add(g);todo.append(g)
 assert len(group)<=60
assert len(group)==60
print(json.dumps({'group_order':60,'status':'PASS','n':17,'k':5,'distance':5,'stabilizer_rank':12,
 'stabilizer_weight_enumerator':swe,'weight_at_most_two_distinct_syndromes':1276,
 'weight_five_logical_witness':witness,'seed':seed,'orbit_size':30,'rows':rows},indent=2))
