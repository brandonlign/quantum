import json,pathlib,itertools,collections,math,hashlib
from fractions import Fraction as F
P=pathlib.Path(__file__).parent;S=json.load(open(P/'generators.json'))['rows'];n=17;m=(1<<n)-1
sp=lambda a,b:(((a&m)&(b>>n)).bit_count()+((a>>n)&(b&m)).bit_count())%2
wt=lambda a:((a&m)|(a>>n)).bit_count()
def need(x,s):
 if not x:raise RuntimeError(s)
need(len(S)==12 and all(0<s<1<<34 for s in S),'rows');need(all(not sp(a,b) for a in S for b in S),'isotropy')
W=[0]
for s in S:W +=[w^s for w in W]
need(len(set(W))==4096,'rank');SS=set(W);A=[0]*18
for w in W:A[wt(w)]+=1
K=lambda i,j:sum((-1)**a*3**(j-a)*math.comb(i,a)*math.comb(n-i,j-a) for a in range(j+1) if a<=i and j-a<=n-i)
Dmw=[]
for j in range(18):
 v=sum(A[i]*K(i,j) for i in range(18));need(v%4096==0,'MacWilliams integrality');Dmw.append(v//4096)
# Independently enumerate the complete normalizer kernel by binary elimination.
B={}
for s in S:
 x=(s>>17)|((s&m)<<17)
 while x:
  j=x.bit_length()-1
  if j in B:x^=B[j]
  else:B[j]=x;break
N=[]
for j in range(34):
 if j in B:continue
 x=1<<j
 for i in sorted(B):
  if (B[i]&x).bit_count()%2:x^=1<<i
 N.append(x)
need(len(N)==22 and all(not sp(v,s) for v in N for s in S),'normalizer basis')
D=[0]*18;D[0]=1;v=0
for i in range(1,1<<22):v^=N[(i&-i).bit_length()-1];D[wt(v)]+=1
need(D==Dmw,'full normalizer / MacWilliams agreement')
# Direct support/color enumeration of every physical error through weight five.
sy=[]
for i in range(17):
 x=sum(sp(1<<i,s)<<j for j,s in enumerate(S));z=sum(sp(1<<(17+i),s)<<j for j,s in enumerate(S));sy.append([x,z,x^z])
low=[];tested=[];witness=None
for w in range(1,6):
 zero=num=0
 for supp in itertools.combinations(range(17),w):
  for colors in itertools.product(range(3),repeat=w):
   h=0;num+=1
   for i,a in zip(supp,colors):h^=sy[i][a]
   if not h:
    zero+=1
    if witness is None:
     e=0
     for i,a in zip(supp,colors):e^=[1<<i,1<<(17+i),(1<<i)|(1<<(17+i))][a]
     witness=e
 low.append(zero);tested.append(num);print('physical weight',w,'tested',num,'zero',zero,flush=True)
need(low==[0,0,0,0,420] and D[1:6]==low,'distance');need(min(i for i in range(1,18) if A[i])==8,'stabilizer minimum');need(witness not in SS,'logical witness')
def g(s):
 v=0
 for k in range(15):
  j=3*(k//3)+(k+1)%3;v|=((s>>k)&1)<<j;v|=((s>>(17+k))&1)<<(17+j)
 for k in [15,16]:
  x=s>>k&1;z=s>>(17+k)&1;v|=z<<k;v|=(x^z)<<(17+k)
 return v
need(all(g(s) in SS for s in S),'C3 symmetry')
# Exact upper bound: distance six would impose five MacWilliams equalities.
y=[F(-1001,201),F(25,25728),F(101,437376),F(197,1749504),F(19,583168),F(5,583168),F(1,583168)]
coeff=[]
for i in range(18):coeff.append(y[0]*int(i==0)+y[1]+sum(y[j+1]*(K(i,j)-4096*int(i==j)) for j in range(1,6)))
need(min(coeff)>=0 and y[0]+4096*y[1]==-1,'exact upper-bound contradiction')
# Non-LC-GF4-linear witness: a stabilizer and normalizer word have only equal-or-zero local labels and odd overlap.
s=S[7];d=8065646656;need(all(not sp(d,a) for a in S),'nonlinearity witness belongs to normalizer');overlap=0
for i in range(17):
 a=((s>>i)&1)|(((s>>(17+i))&1)<<1);b=((d>>i)&1)|(((d>>(17+i))&1)<<1);need(not a or not b or a==b,'local commuting witness');overlap+=bool(a and b)
need(overlap==3,'odd overlap')
names='IXZY';strings=[''.join(names[((s>>i)&1)|(((s>>(17+i))&1)<<1)] for i in range(17)) for s in S];(P/'generators.txt').write_text('\n'.join('+'+x for x in strings)+'\n')
out={'status':'PASS','rank':12,'isotropic':True,'pure_distance':5,'stabilizer_minimum':8,'stabilizer_enumerator':A,'normalizer_enumerator':D,'physical_errors_tested':tested,'physical_zero_syndromes':low,'weight5_logical_witness':witness,'C3_invariant':True,'distance_optimal_among_binary_stabilizers':True,'not_LC_equivalent_to_GF4_linear':True,'nonlinearity_witness':[s,d],'upper_bound_coefficients':list(map(str,coeff))};(P/'VERIFICATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
