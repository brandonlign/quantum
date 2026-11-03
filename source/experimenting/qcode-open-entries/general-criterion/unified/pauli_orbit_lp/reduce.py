"""Eliminate exactly fixed nonnegative variables, retaining certificate lifts."""
from fractions import Fraction
import numpy as np

def reduced(E,b,G,identity=0):
 fixed={};unit_rows={}
 for i,row in enumerate(E):
  nz=np.flatnonzero(row)
  if len(nz)==1:
   j=int(nz[0]);v=Fraction(int(b[i]),int(row[j]));assert v>=0
   if j in fixed:assert fixed[j]==v
   fixed[j]=v;unit_rows[j]=i
 assert fixed[identity]==1
 fixed.pop(identity);active=[j for j in range(E.shape[1]) if j not in fixed];ii=active.index(identity)
 assert all(v.denominator==1 for v in fixed.values())
 def transform(M):
  out=M[:,active].copy()
  for j,v in fixed.items():out[:,ii]+=int(v)*M[:,j]
  return out
 er=transform(E);gr=transform(G);keep=[]
 for i,row in enumerate(er):
  if np.any(row):keep.append(i)
  else:assert b[i]==0
 return er[keep],b[keep],gr,dict(active=active,fixed=fixed,unit_rows=unit_rows,kept_equations=keep,identity=identity)

def lift_certificate(nu,lam,E,b,G,meta):
 out=[Fraction(0)]*len(E)
 for i,v in zip(meta['kept_equations'],nu):out[i]+=v
 coeff=[sum(v*int(E[i,j]) for i,v in enumerate(out) if v)+sum(v*int(G[i,j]) for i,v in enumerate(lam) if v) for j in range(E.shape[1])]
 idrow=meta['unit_rows'][meta['identity']];idcoef=int(E[idrow,meta['identity']]);assert idcoef==1
 for j,value in meta['fixed'].items():
  alpha=-coeff[j];u=meta['unit_rows'][j];a=int(E[u,j]);out[u]+=alpha/a;out[idrow]-=alpha*value
 return out,lam
