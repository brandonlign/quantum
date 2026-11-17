"""Self-contained exact verification of derivative code_b_style: pure [[18,3,6]].
Pauli vector convention: x | (z << 18). No third-party packages or files.
"""
from collections import Counter
from math import comb
import json
if not __debug__:raise RuntimeError('Run without -O: assertions are proof checks')
S=[30539046912, 59444428930, 57231732741, 7398514697, 14877519889, 27939967136, 16496890048, 63874044160, 17245201025, 66477130752, 15987673088, 25274630273, 51785547776, 16911630336, 66846720]
mask=(1<<18)-1
def wt(v):return ((v&mask)|(v>>18)).bit_count()
def pair(a,b):return (((a&mask)&(b>>18))^((b&mask)&(a>>18))).bit_count()%2
def rref(rows,width):
 rows=rows[:];piv=[];i=0
 for j in range(width):
  k=next((k for k in range(i,len(rows)) if rows[k]>>j&1),None)
  if k is None:continue
  rows[i],rows[k]=rows[k],rows[i]
  for k in range(len(rows)):
   if k!=i and rows[k]>>j&1:rows[k]^=rows[i]
  piv.append(j);i+=1
 return rows[:i],piv
def words(rows):
 x=0;yield x
 for k in range(1,1<<len(rows)):
  x^=rows[(k&-k).bit_length()-1];yield x
assert len(S)==15 and all(0<v<1<<36 for v in S)
assert len(rref(S,36)[0])==15 and all(pair(a,b)==0 for a in S for b in S)
equations=[(s>>18)|((s&mask)<<18) for s in S];rr,piv=rref(equations,36)
N=[(1<<j)^sum(((r>>j)&1)<<p for r,p in zip(rr,piv)) for j in range(36) if j not in piv]
assert len(N)==len(rref(N,36)[0])==21 and all(pair(s,n)==0 for s in S for n in N)
stab=set(words(S));assert len(stab)==1<<15
A=Counter(wt(s) for s in stab);B=Counter();L=Counter();witness=None
for v in words(N):
 w=wt(v);B[w]+=1
 if v not in stab:
  L[w]+=1
  if w==6 and witness is None:witness=v
assert sum(B.values())==1<<21 and min(w for w in B if w)==6 and min(L)==6
assert all(B[w]==A[w]+L[w] for w in range(19))
# Independent q=4 MacWilliams transform of the stabilizer enumerator.
for j in range(19):
 value=sum(A[i]*sum((-1)**a*3**(j-a)*comb(i,a)*comb(18-i,j-a) for a in range(max(0,j-(18-i)),min(i,j)+1)) for i in range(19))
 assert value==(1<<15)*B[j]
assert min(w for w in A if w)==7 and A[6]==0 and A[7]==16 and B[6]==564 and L[6]==564
explicit=34359943556
assert wt(explicit)==6 and explicit not in stab and all(pair(explicit,s)==0 for s in S)
print(json.dumps({'status':'PASS','parameters':'pure [[18,3,6]]','rank_S':15,'rank_normalizer':21,'normalizer_words_enumerated':sum(B.values()),'minimum_normalizer_weight':6,'logical_distance':6,'explicit_weight6_logical':explicit,'stabilizer_enumerator':dict(sorted(A.items())),'normalizer_enumerator':dict(sorted(B.items())),'logical_enumerator':dict(sorted(L.items())),'all_19_MacWilliams_coefficients_match':True},indent=2))
