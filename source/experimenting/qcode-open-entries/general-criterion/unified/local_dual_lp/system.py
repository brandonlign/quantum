"""Necessary letter-refined split LP with both local centralizers fixed."""
import sys
from pathlib import Path
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];EXPERIMENT=ROOT.parents[1]
sys.path[:0]=[str(EXPERIMENT/'pure'),str(EXPERIMENT/'split'),str(ROOT)]
import refined_split as rs
from model import span,weight

def letter_hist(words,m):
 return Counter(((w>>m)&~w&((1<<m)-1)).bit_count()*100+(w&((1<<m)-1)).bit_count() for w in words)

def build(n,k,d,m,Hhist,Chist,Nhist,light_budget=None,weight8_min=None):
 E,b,G,V,*_=rs.build(m,d,n=n,k=k,d=d,pure=True);M=1<<(n-k);f=n-m;VI={v:i for i,v in enumerate(V)}
 bk={(z,t):rs.pmul(rs.pmul(rs.ppow(rs.PI,m-z-t),rs.ppow(rs.PZ,z)),rs.ppow(rs.PX,t)) for t in range(m+1) for z in range(m-t+1)}
 def members(key):
  z,t=key;return sorted({(z,t),(m-z-t,t)})
 def dual(zc,tc,bc):
  return [sum(bk[(z,t)].get((zc,tc),0)*rs.kraw(f,bc,w) for z,t in members(key)) for key,w in V]
 rows=[];rhs=[];labels=[]
 def add(row,value,label):rows.append(row);rhs.append(value);labels.append(label)
 # S intersect prefix is exactly <Z^m>.
 for key,w in V:
  if w:continue
  row=[0]*len(V);row[VI[(key,w)]]=1
  add(row,int(key==(0,0)),('local_S',key))
 # S intersect suffix is exactly the fixed N.
 for w in range(f+1):
  row=[0]*len(V);row[VI[((0,0),w)]]=1;add(row,Nhist.get(w,0),('suffix_S',w))
 # D intersect prefix is exactly the local logical code H (doubled by R).
 for t in range(m+1):
  for z in range(m-t+1):add(dual(z,t,0),M*Hhist.get(z*100+t,0),('local_D',z,t))
 for w in range(f+1):add(dual(0,0,w),M*Chist.get(w,0),('suffix_D',w))
 extra=[];identity=[0]*len(V);identity[VI[((0,0),0)]]=1
 def aw(w):return [sum(z+t+b==w for z,t in members(key)) for key,b in V]
 if light_budget is not None:extra.append([x+y-light_budget*i for x,y,i in zip(aw(6),aw(7),identity)])
 if weight8_min is not None:extra.append([weight8_min*i-x for i,x in zip(identity,aw(8))])
 return np.vstack([E,np.array(rows,dtype=object)]),np.concatenate([b,np.array(rhs,dtype=object)]),np.vstack([G,np.array(extra,dtype=object)]) if extra else G,V,labels
