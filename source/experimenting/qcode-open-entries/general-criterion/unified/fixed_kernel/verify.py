"""Independent complete fixed-form trees and physical fixed-kernel Y maps."""
import itertools,json,sys,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,gram,lift_coordinate,reduced_basis,span,weight,coordinates,symp
from fixed_suffix_search import quotient_input
from fixed_kernel.common import standard,prefix_data,payload,kernel
from fixed_kernel.solve import solve
from mixed_frames.fixed_z import solve as reference
from mixed_frames.stage3.common import atomic_json,sha

def main():
 fixture=json.loads((HERE/'control.json').read_text());H=fixture['H'];model=ExtensionModel(8,2,3,4,[15<<4],H);data=quotient_input(model,[15<<4]);prefix=prefix_data(H,4);reps=[v for v,d in enumerate(data['distances']) if d>=2]
 binary=HERE/'z_enum_fixed';cmd=[str(binary),'--gamma','--avoid-radical','--first-orbits','--fixed-form'];p=subprocess.run(cmd,input=payload(model,data,reps,prefix),text=True,capture_output=True,check=True);stats=json.loads(p.stderr);L=[tuple(map(int,l.split()[1:])) for l in p.stdout.splitlines()];assert stats['status']=='COMPLETE' and stats['leaves']==len(L)==len(set(L))
 G=data['Gamma'];inverse={lift_coordinate(c,G):c for c in range(16)};expected=[]
 for xs in itertools.product(reps,repeat=3):
  if len(reduced_basis(xs))!=3:continue
  if any(((inverse[xs[i]]&xs[j]).bit_count()%2)!=((prefix['Gram'][i]>>j)&1) for i in range(3) for j in range(3)):continue
  if any(data['distances'][lift_coordinate(a,xs)]<3-min(a.bit_count(),4-a.bit_count()) for a in range(1,8)):continue
  expected.append(xs)
 assert set(L)==set(expected)
 frame=standard(4);ours=set();theirs=set();positive=0
 for xs in L:
  full=[lift_coordinate(x,data['units']) for x in xs];a=solve(model,H,full,find_all=True);b=reference(model,frame,full,find_all=True)
  assert a['status']==b['status']=='COMPLETE'
  mine={tuple(s['standard_columns']) for s in a['solutions']}
  ref=set()
  for s in b['solutions']:
   columns=s['adapted_columns']
   if all(model.map_value(coordinates(h,frame.R+model.E)>>1,frame.model_columns(model,columns))==0 for h in H):ref.add(tuple(columns))
  assert mine==ref
  ours.update(mine);theirs.update(ref)
  for s in a['solutions']:
   S,D=model.construct(s['model_columns']);assert min(weight(w,8) for w in span(D) if w)>=3
   from gauge_fix import supported_subspace
   assert supported_subspace(D,8,range(4))==reduced_basis([__import__('model').embed(h,0,4,4) for h in H])
   positive+=1
 assert ours==theirs and positive
 group=HERE/'control_group.txt';group.write_text(f'{len(prefix["group"])} 3\n'+'\n'.join(' '.join(map(str,g)) for g in prefix['group'])+'\n')
 p=subprocess.run(cmd+['--prefix-group',str(group)],input=payload(model,data,reps,prefix),text=True,capture_output=True,check=True);s=json.loads(p.stderr);K=[tuple(map(int,l.split()[1:])) for l in p.stdout.splitlines()];assert s['status']=='COMPLETE'
 canonical={min(tuple(lift_coordinate(u,xs) for u in g) for g in prefix['group']) for xs in L};assert set(K)==canonical
 for limit in (0,1):
  p=subprocess.run(cmd+['--max-nodes',str(limit)],input=payload(model,data,reps,prefix),text=True,capture_output=True,check=True);assert json.loads(p.stderr)['status']=='LIMIT'
 assert solve(model,H,[lift_coordinate(x,data['units']) for x in L[0]],node_limit=0)['status']=='LIMIT'
 report=dict(status='PASS',parameters=[8,2,3],complete_Z_leaves=len(L),prefix_representatives=len(K),physical_prefix_group=prefix['physical_order'],Z_group=len(prefix['group']),full_physical_maps=positive,reference_map_set_matches=True,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'z_enum_fixed.c',binary,HERE/'common.py',HERE/'solve.py',ROOT/'model.py',ROOT/'mixed_frames/fixed_z.py',ROOT/'fixed_suffix_search.py']})
 atomic_json(HERE/'verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
