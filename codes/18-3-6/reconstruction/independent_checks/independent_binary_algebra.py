import hashlib,json,random,pathlib
def sp(a,b,n):
 m=(1<<n)-1
 return (((a&m)&(b>>n)).bit_count()+((b&m)&(a>>n)).bit_count())&1
def wt(a,n):return ((a&((1<<n)-1))|(a>>n)).bit_count()
def span(rows):
 out=[0]
 for r in rows:out +=[v^r for v in out]
 return out
def basis(rows):
 d={}
 for v in rows:
  while v:
   k=v.bit_length()-1
   if k in d:v^=d[k]
   else:d[k]=v;break
 return d
# Independent Gauss-Jordan affine solver, bit equations dot(x,row)=rhs.
def solve(rows,rhs,n):
 rows=[r|(s<<n) for r,s in zip(rows,rhs)];piv=[];i=0
 for j in range(n):
  k=next((k for k in range(i,len(rows)) if rows[k]>>j&1),None)
  if k is None:continue
  rows[i],rows[k]=rows[k],rows[i]
  for k in range(len(rows)):
   if k!=i and rows[k]>>j&1:rows[k]^=rows[i]
  piv.append(j);i+=1
 assert all((r&((1<<n)-1)) or not(r>>n) for r in rows)
 x=sum((rows[k]>>n)<<j for k,j in enumerate(piv))
 return x,rows[:i],piv
def perp(rows,n):
 eq=[(r>>n)|((r&((1<<n)-1))<<n) for r in rows]
 _,rr,piv=solve(eq,[0]*len(eq),2*n)
 return [(1<<j)^sum(((r>>j)&1)<<k for r,k in zip(rr,piv)) for j in range(2*n) if j not in piv]
def emb(v,n,off):return ((v&((1<<n)-1))<<off)|((v>>n)<<(18+off))
def sp8(a,b):return sum(((a>>(2*i)&1)*(b>>(2*i+1)&1))^((b>>(2*i)&1)*(a>>(2*i+1)&1)) for i in range(4))&1
