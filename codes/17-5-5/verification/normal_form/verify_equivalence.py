import json,pathlib,collections
P=pathlib.Path(__file__).resolve().parent;O=P;OLD=P.parent
if not __debug__:raise RuntimeError('Run without Python optimization flags')
f=json.load(open(P/'NORMAL_FORM.json'));c=json.load(open(OLD/'a5_construction.json'));orig=json.load(open(OLD/'generators.json'))['rows'];n=17;mask=(1<<n)-1
I=(0,1,2,3)
def sp(a,b):return (((a&mask)&(b>>n)).bit_count()+((b&mask)&(a>>n)).bit_count())&1
def weight(a):return ((a&mask)|(a>>n)).bit_count()
def span(rows):
 s={0}
 for v in rows:s|={w^v for w in s}
 return s
def rref(rows):
 rows=list(rows);p=0
 for j in range(34):
  k=next((k for k in range(p,len(rows)) if rows[k]>>j&1),None)
  if k is None:continue
  rows[k],rows[p]=rows[p],rows[k]
  for k in range(len(rows)):
   if k!=p and rows[k]>>j&1:rows[k]^=rows[p]
  p+=1
 return rows[:p]
def monomial(g):
 p=tuple(g['permutation']);m=tuple(tuple(v) for v in g['local_label_maps'])
 assert sorted(p)==list(range(n)) and len(m)==n
 assert all(v[0]==0 and sorted(v)==[0,1,2,3] and all(v[x^y]==v[x]^v[y] for x in range(4) for y in range(4)) for v in m)
 return p,m
def act(g,v):
 p,m=g;out=0
 for i in range(n):
  x=((v>>i)&1)|(((v>>(n+i))&1)<<1);y=m[i][x];out|=(y&1)<<p[i]|(y>>1)<<(n+p[i])
 return out
def multiply(g,h):
 # Apply h first, then g; local maps are attached to source sites.
 pg,mg=g;ph,mh=h
 return tuple(pg[ph[i]] for i in range(n)),tuple(tuple(mg[ph[i]][mh[i][x]] for x in range(4)) for i in range(n))
identity=tuple(range(n)),tuple(I for _ in range(n))
def order(g):
 h=identity
 for k in range(1,121):
  h=multiply(g,h)
  if h==identity:return k
 raise AssertionError('order > 120')
new=list(map(monomial,f['generators']));old=[monomial(c['generator_'+x]) for x in 'ab'];a,b=new
# Construct new -> old as a binary matrix, then invert by linear elimination.
T=monomial({'permutation':f['coordinate_order'],'local_label_maps':f['local_new_to_old']})
cols=[act(T,1<<i) for i in range(34)];assert all(sp(cols[i],cols[j])==sp(1<<i,1<<j) for i in range(34) for j in range(34));echelon={}
for j,v in enumerate(cols):
 lab=1<<j
 while v:
  p=v.bit_length()-1
  if p in echelon:v^=echelon[p][0];lab^=echelon[p][1]
  else:echelon[p]=(v,lab);break
 assert v
assert len(echelon)==34
def old_to_new(v):
 out=0
 while v:
  p=v.bit_length()-1;x,lab=echelon[p];v^=x;out^=lab
 return out
for v in [1<<i for i in range(34)]:
 assert old_to_new(act(T,v))==v
 for ng,og in zip(new,old):assert act(T,act(ng,v))==act(og,act(T,v))
assert f['stabilizer_rows']==[old_to_new(v) for v in orig]
assert f['seed']==old_to_new(c['seed'])
word=''.join('IXZY'[((f['seed']>>i)&1)|(((f['seed']>>(n+i))&1)<<1)] for i in range(n));assert word==f['seed_paulis']=='IIZXYIYXIZIIZYIII'
S=span(f['stabilizer_rows']);assert len(S)==4096 and all(not sp(u,v) for u in f['stabilizer_rows'] for v in f['stabilizer_rows'])
assert S=={old_to_new(v) for v in span(orig)}
for g in new:assert all(act(g,v) in S for v in f['stabilizer_rows'])
assert all(v==I for v in a[1]);twists=[i for i,v in enumerate(b[1]) if v!=I];assert twists==[14,16] and all(b[0][i]==i for i in twists)
assert b[1][14]==(0,2,3,1) and b[1][16]==(0,3,1,2)
assert all(b[1][14][b[1][16][x]]==x for x in range(4))
assert [order(a),order(b),order(multiply(a,b))]==[2,3,5]
group={identity};queue=[identity]
for h in queue:
 for g in new:
  gh=multiply(g,h)
  if gh not in group:group.add(gh);queue.append(gh)
 assert len(group)<=60
assert len(group)==60
orbits=[];unused=set(range(n))
while unused:
 i=min(unused);orb={g[0][i] for g in group};orbits.append(sorted(orb));unused-=orb
assert orbits==[list(range(12)),list(range(12,17))]
restrictions={tuple(g[0][i]-12 for i in range(12,17)) for g in group};assert len(restrictions)==60
assert all(sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0 for p in restrictions)
orbit={act(g,f['seed']) for g in group};assert len(orbit)==30 and all(weight(v)==8 for v in orbit) and span(orbit)==S
# Independently close the old orbit, then compare its transformed members.
oldorbit={c['seed']};todo=[c['seed']]
for v in todo:
 for g in old:
  w=act(g,v)
  if w not in oldorbit:oldorbit.add(w);todo.append(w)
assert {old_to_new(v) for v in oldorbit}==orbit
A=[0]*18
for v in S:A[weight(v)]+=1
assert min(weight(v) for v in S if v)==8
out={'status':'PASS','conjugacy_basis_checks':68,'original_stabilizer_rows_matched':12,'all_transformed_codewords_matched':4096,'rank':12,'isotropic':True,'seed_paulis':word,'group_order':60,'generator_and_product_orders':[2,3,5],'group_identification':'faithful even-permutation action on five sites, hence A5','site_orbits':orbits,'orbit_size':30,'orbit_weights':[8],'orbit_span_equals_transformed_code':True,'nontrivial_local_maps_b':twists,'stabilizer_enumerator':A,'scope':'Equivalent binary monomial presentation only; no new code/family/full automorphism claim or phase-consistent unitary A5 claim.'}
(O/'NORMAL_FORM_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');(O/'rows.txt').write_text('\n'.join(map(str,f['stabilizer_rows']))+'\n');(O/'orbit.txt').write_text('\n'.join(map(str,sorted(orbit)))+'\n')
print(json.dumps(out,indent=2))
