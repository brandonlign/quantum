"""Full local Pauli-orbit split relaxation; character sums computed exactly."""
from pathlib import Path
import sys
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'m8r1'),str(ROOT.parents[1]/'split')]
from model import coordinates,dual,span,weight,symp,reduced_basis,lift_coordinate
from mixed_frames.frames import xor_table
from suffix_symmetries import automorphisms
from system import kraw

def complement(space,sub):
 seen=list(sub);out=[]
 for w in space:
  if len(reduced_basis(seen+[w]))>len(seen):seen.append(w);out.append(w)
 return out

def orbit_partition(columns,width):
 N=1<<width;parent=list(range(N))
 def find(v):
  while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
  return v
 for cols in columns:
  for a,b in enumerate(xor_table(cols)):
   x,y=find(a),find(b)
   if x!=y:parent[max(x,y)]=min(x,y)
 groups={}
 for v in range(N):groups.setdefault(find(v),[]).append(v)
 return list(groups.values())

def fwht(values):
 values=values.copy();h=1
 while h<len(values):
  blocks=values.reshape(-1,2*h);a=blocks[:,:h].copy();b=blocks[:,h:].copy();blocks[:,:h]=a+b;blocks[:,h:]=a-b;h*=2
 return values

def prepare(H,m):
 H=reduced_basis(H);p=((1<<m)-1)<<m
 assert m%2==0 and p in span(H) and all(symp(p,h,m)==0 for h in H)
 A=complement(dual(H,m),[p]);B=complement(dual([p],m),H);width=len(A);assert len(B)==width
 labels=[sum(symp(a,b,m)<<i for i,a in enumerate(A)) for b in B]
 assert len(reduced_basis(labels))==width
 B=[lift_coordinate(coordinates(1<<i,labels),B) for i in range(width)]
 assert all(symp(a,b,m)==(i==j) for i,a in enumerate(A) for j,b in enumerate(B))
 maps,order=automorphisms(H,m,[p]);Ac=[];Bc=[]
 for apply in maps:
  Ac.append([coordinates(apply(a),[p]+A)>>1 for a in A])
  Bc.append([sum(symp(a,apply(b),m)<<i for i,a in enumerate(A)) for b in B])
 AO=orbit_partition(Ac,width);BO=orbit_partition(Bc,width);aw=xor_table(A);bw=xor_table(B);hs=span(H)
 char=np.zeros((len(BO),len(AO)),dtype=np.int64)
 for j,orbit in enumerate(AO):
  values=np.zeros(1<<width,dtype=np.int64);values[orbit]=1;values=fwht(values)
  for i,ob in enumerate(BO):
   assert np.all(values[ob]==values[ob[0]]);char[i,j]=values[ob[0]]
 weights=[];sign=[]
 for orbit in AO:
  pairs={tuple(sorted((weight(aw[a],m),weight(aw[a]^p,m)))) for a in orbit}
  parity={weight(aw[a],m)%2 for a in orbit};assert len(pairs)==len(parity)==1
  pair=next(iter(pairs));assert pair[0]%2==pair[1]%2;weights.append(pair);sign.append((-1)**next(iter(parity)))
 costs=[]
 for orbit in BO:
  minima={min(weight(bw[b]^h,m) for h in hs) for b in orbit};assert len(minima)==1;costs.append(next(iter(minima)))
 blookup={v:i for i,o in enumerate(BO) for v in o};incl=[]
 for orbit in AO:
  images={blookup[sum(symp(a,aw[v],m)<<j for j,a in enumerate(A))] for v in orbit};assert len(images)==1;incl.append(next(iter(images)))
 return dict(m=m,H=H,p=p,A=A,B=B,AO=AO,BO=BO,char=char,weights=weights,sign=sign,costs=costs,incl=incl,physical_order=order)

def build(data,n,k,d,C_hist,N_hist,light_budget=None,weight8_min=None):
 m=data['m'];f=n-m;ell=n-k-(2*m-len(data['H']));assert ell>=0 and sum(N_hist.values())==1<<ell
 AO=data['AO'];BO=data['BO'];na=len(AO);nb=len(BO);F=f+1;den=1<<(n-k-1);cols=na*F
 K=np.array([[kraw(f,j,b) for b in range(F)] for j in range(F)],dtype=np.int64)
 T=np.kron(data['char'],K)
 eq=[];rhs=[];le=[];labels=[]
 def E(row,value,label):eq.append(row);rhs.append(value);labels.append(label)
 def unit(a,b):
  row=np.zeros(cols,dtype=np.int64);row[a*F+b]=1;return row
 remaining_light=None if light_budget is None else light_budget-N_hist.get(6,0)-N_hist.get(7,0)
 if remaining_light is not None:assert remaining_light>=0
 for a,o in enumerate(AO):
  row=np.zeros(cols,dtype=np.int64);row[a*F:(a+1)*F]=1;E(row,1<<ell,('prefix_marginal',a))
  for b in range(F):
   if a==0:E(unit(a,b),N_hist.get(b,0),('suffix_S',b))
   elif (b==0 or min(data['weights'][a])+b<d or
         (remaining_light is not None and sum(v+b in (6,7) for v in data['weights'][a])>remaining_light)):
    E(unit(a,b),0,('local_pure_or_integer_light_S',a,b))
 for c,o in enumerate(BO):
  for j in range(F):
   row=T[c*F+j]
   if c==0:E(row,den*C_hist.get(j,0),('suffix_D',j))
   elif j==0 or data['costs'][c]+j<d:E(row,0,('local_or_pure_D',c,j))
   le.append(-row)
 signs=np.repeat(np.array(data['sign'],dtype=np.int64),F)*np.tile(np.array([(-1)**b for b in range(F)],dtype=np.int64),na)
 le.extend(-T*signs)
 for a in range(na):
  c=data['incl'][a]
  for b in range(F):le.append(den*unit(a,b)-T[c*F+b])
 def aw(w):return np.array([len(AO[a])*sum(v+b==w for v in data['weights'][a]) for a in range(na) for b in range(F)],dtype=np.int64)
 identity=unit(0,0)
 if light_budget is not None:le.append(aw(6)+aw(7)-light_budget*identity)
 if weight8_min is not None:le.append(weight8_min*identity-aw(8))
 return np.array(eq,dtype=np.int64),np.array(rhs,dtype=np.int64),np.array(le,dtype=np.int64),labels
