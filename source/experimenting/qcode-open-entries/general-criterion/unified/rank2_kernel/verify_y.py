"""Independent full-coordinate linear-system oracle for fixed-Z Y maps."""
import json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,lift_coordinate,reduced_basis,span,weight
from rank2_kernel.common import local_R,frame_for,kernel,quotient_input
from rank2_kernel.solve import solve
from mixed_frames.fixed_z import solve as reference
from mixed_frames.stage3.common import atomic_json,sha

def affine(rows,width):
 piv={}
 for v,rhs in rows:
  for j in sorted(piv,reverse=True):
   if v>>j&1:v^=piv[j][0];rhs^=piv[j][1]
  if v:piv[v.bit_length()-1]=(v,rhs)
  elif rhs:return None
 def complete(v,homogeneous):
  for j,(row,rhs) in sorted(piv.items()):
   if ((row&v).bit_count()%2)^(0 if homogeneous else rhs):v^=1<<j
  return v
 base=complete(0,False);generators=[complete(1<<i,True) for i in range(width) if i not in piv]
 assert all((mask&base).bit_count()%2==rhs for mask,rhs in rows)
 assert all(all((mask&v).bit_count()%2==0 for mask,rhs in rows) for v in generators)
 return base,generators

def oracle(model,H,frame,X,N):
 p=frame.p;b=model.b;rows=[]
 for i in range(b):
  for j in range(i):
   mask=sum(((x>>i)&1)<<(a*b+j) | ((x>>j)&1)<<(a*b+i) for a,x in enumerate(X))
   rows.append((mask,(model.Gamma[i]>>j)&1))
 for h in kernel(H,frame.R,model.m):
  e=coordinates(h,frame.R+frame.E)>>2;z=e&((1<<p)-1);x=e>>p;rhs=lift_coordinate(z,X)
  for j in range(b):rows.append((sum(((x>>a)&1)<<(a*b+j) for a in range(p)),(rhs>>j)&1))
 for n in N:
  coeff=coordinates(n,model.P)
  for a in range(p):rows.append((coeff<<(a*b),0))
 answer=affine(rows,p*b);assert answer is not None;base,G=answer;assert len(G)<=10
 result=set()
 for delta in span(G):
  v=base^delta;Y=[(v>>(a*b))&((1<<b)-1) for a in range(p)];columns=frame.model_columns(model,X+Y)
  if model.criterion(columns):result.add(tuple(columns))
 return result,len(G)

def main():
 start=time.monotonic();out=[]
 for r in json.loads((HERE/'controls.json').read_text())['records']:
  n,k,d=r['parameters'];model=ExtensionModel(n,k,d,r['m'],local_R(r['m']),r['C']);frame=frame_for(r['H'],r['m']);adapted=frame.adapted_columns(model,r['model_columns']);X=adapted[:frame.p]
  expected,free=oracle(model,r['H'],frame,X,r['N']);assert expected
  answer=solve(model,r['H'],frame,X,find_all=True);assert answer['status']=='COMPLETE'
  actual={tuple(s['model_columns']) for s in answer['solutions']};assert actual==expected
  if n<=9:
   old=reference(model,frame,X,find_all=True);assert old['status']=='COMPLETE'
   oldset={tuple(s['model_columns']) for s in old['solutions'] if len(reduced_basis(s['model_columns']))==2*frame.p-(len(r['H'])-2) and all(model.map_value(coordinates(h,frame.R+model.E)>>2,s['model_columns'])==0 for h in r['H'])}
   assert oldset==expected
  assert solve(model,r['H'],frame,X,node_limit=0)['status']=='LIMIT'
  assert solve(model,r['H'],frame,X,seconds=0)['status']=='LIMIT'
  out.append(dict(parameters=r['parameters'],q=len(r['H'])-2,p=frame.p,free_linear_bits=free,complete_physical_maps=len(actual),CSP_nodes=answer['nodes'],reference_compared=n<=9))
 report=dict(status='PASS',records=out,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'solve.py',HERE/'common.py',HERE/'controls.json',ROOT/'model.py',ROOT/'two_site.py',ROOT/'mixed_frames/fixed_z.py']})
 atomic_json(HERE/'y_verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
