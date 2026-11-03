"""Exact fixed-Z CSP with all local-kernel equations, no phase quotient."""
import time
import numpy as np
from model import coordinates,lift_coordinate,reduced_basis
from mixed_frames.frames import xor_table
from two_site import normal_form
from fixed_kernel.common import standard,kernel

class Limit(Exception):pass

def solve(model,H,xs,*,find_all=False,node_limit=None,seconds=None):
 start=time.monotonic();frame=standard(model.m);p=frame.p;Q=1<<p;nodes=0;solutions=[];reason=None
 assert len(xs)==p and len(reduced_basis(xs))==p
 def check():
  if seconds is not None and time.monotonic()-start>=seconds:raise Limit
 def result(status,complete):return dict(status=status,complete=complete,fixed_z_excluded=complete and not solutions,nodes=nodes,reason=reason,solutions=solutions,seconds=time.monotonic()-start)
 try:
  check()
  if node_limit==0:raise Limit
  nf=normal_form(xs,model.Gamma);G=nf['Gamma'];forward=np.asarray(nf['forward'],dtype=np.int64)
  if any(G[i]>>p for i in range(p,model.b)):
   reason='annihilator_Gram';return result('COMPLETE',True)
  zsyn=np.asarray(xor_table(xs),dtype=np.int64);offset=xor_table(nf['offsets']);dist=np.asarray(model.distances,dtype=np.int8)
  if any(frame.costs[z]+model.distances[int(zsyn[z])]<model.d for z in range(1,Q)):
   reason='pure_Z_weight';return result('COMPLETE',True)
  equations=[]
  for h in kernel(H,model.m):
   c=coordinates(h,frame.R+frame.E)>>1;z=c&(Q-1);mask=c>>p;rhs=z^lift_coordinate(mask,nf['offsets'])
   if rhs>=Q:
    reason='kernel_projection';return result('COMPLETE',True)
   equations.append((mask,rhs))
  allowed=[0]*Q;u=np.arange(Q,dtype=np.int64)
  for sector in range(1,Q):
   check();cost=np.asarray(frame.costs[sector*Q:(sector+1)*Q]);requirements=np.flatnonzero(cost<model.d);ys=forward[u^offset[sector]]
   ok=np.all(dist[ys[:,None]^zsyn[requirements][None,:]]>=model.d-cost[requirements][None,:],axis=1)
   allowed[sector]=sum(1<<int(v) for v in np.flatnonzero(ok))
   if not allowed[sector]:reason='empty_sector_'+str(sector);return result('COMPLETE',True)
  domains=[allowed[1<<i] for i in range(p)];values=[None]*p
  parity_masks=[[sum(1<<v for v in range(Q) if (v>>i)&1==want) for want in range(2)] for i in range(p)]
  def translate(bits,shift):
   out=0
   while bits:
    low=bits&-bits;out|=1<<((low.bit_length()-1)^shift);bits^=low
   return out
  def dfs(dom,assigned):
   nonlocal nodes
   check()
   if node_limit is not None and nodes>=node_limit:raise Limit
   nodes+=1;dom=dom[:]
   for sector in range(1,Q):
    remaining=sector&~assigned;done=sector&assigned;total=0
    while done:
     bit=done&-done;total^=values[bit.bit_length()-1];done^=bit
    if not remaining:
     if not (allowed[sector]>>total)&1:return False
    elif remaining&(remaining-1)==0:
     i=remaining.bit_length()-1;dom[i]&=translate(allowed[sector],total)
     if not dom[i]:return False
   for mask,rhs in equations:
    remaining=mask&~assigned;done=mask&assigned;total=rhs
    while done:
     bit=done&-done;total^=values[bit.bit_length()-1];done^=bit
    if not remaining:
     if total:return False
    elif remaining&(remaining-1)==0:
     i=remaining.bit_length()-1;dom[i]&=1<<total
     if not dom[i]:return False
   if assigned==Q-1:
    ys=[int(forward[values[i]^nf['offsets'][i]]) for i in range(p)];columns=frame.model_columns(model,list(xs)+ys)
    assert model.criterion(columns)
    assert all(model.map_value(coordinates(h,frame.R+model.E)>>1,columns)==0 for h in H)
    solutions.append(dict(model_columns=columns,standard_columns=list(xs)+ys,U=values[:]));return not find_all
   i=min((i for i in range(p) if not assigned>>i&1),key=lambda i:dom[i].bit_count());options=dom[i]
   while options:
    bit=options&-options;options^=bit;v=bit.bit_length()-1;values[i]=v;child=dom[:];ok=True
    for j in range(p):
     if j==i or assigned>>j&1:continue
     want=((G[i]>>j)^(v>>j))&1;child[j]&=parity_masks[i][want]
     if not child[j]:ok=False;break
    if ok and dfs(child,assigned|(1<<i)):return True
   return False
  found=dfs(domains,0)
  return result('FOUND' if found else 'COMPLETE',not found)
 except Limit:reason='search_limit';return result('LIMIT',False)
