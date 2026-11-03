"""Complete small fixed-kernel trees and independently filtered physical maps."""
import itertools,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,dual,gram,lift_coordinate,nullspace,reduced_basis,span,split,symp,weight
from fixed_suffix_search import quotient_input
from gauge_fix import supported_subspace
from fixed_kernel.common import kernel_contexts,payload,standard
from fixed_kernel.solve import solve
from mixed_frames.fixed_z import solve as reference
from mixed_frames.stage3.common import atomic_json,sha
from run_benchmark import z_condition

def native(model,data,context,*,group=None,limit=None,fixed=None):
 reps=[v for v,d in enumerate(data['distances']) if d>=model.d-1]
 cmd=[str(HERE/'z_enum_fixed'),'--gamma','--avoid-radical','--first-orbits','--kernel-data']
 if group is not None:cmd+=['--prefix-group',str(group)]
 if limit is not None:cmd+=['--max-nodes',str(limit)]
 if fixed is not None:cmd+=['--prefix',','.join(map(str,fixed))]
 proc=subprocess.run(cmd,input=payload(model,data,reps,context),capture_output=True,text=True,check=True,timeout=60)
 stats=json.loads(proc.stderr);leaves=[tuple(map(int,l.split()[1:])) for l in proc.stdout.splitlines()]
 assert len(leaves)==len(set(leaves))==stats['leaves']
 return stats,leaves

def oracle(model,data,c):
 p=model.m-1;g=c['g'];Gamma=data['Gamma'];pre={lift_coordinate(a,Gamma):a for a in range(1<<len(Gamma))}
 reps=[v for v,d in enumerate(data['distances']) if d>=model.d-1];out=[]
 for xs in itertools.product(reps,repeat=p):
  if len(reduced_basis(xs))!=p:continue
  if any(((r&xs[i]).bit_count()%2)!=((c['tags'][i]>>j)&1) for i in range(p) for j,r in enumerate(c['radical_coeff'])):continue
  ts=[xs[i]^lift_coordinate(coordinates(c['tags'][i],c['tags'][:g]),xs[:g]) for i in range(p)]
  assert all(t in pre for t in ts)
  if any(((pre[ts[i]]&ts[j]).bit_count()%2)!=((c['Gram'][i]>>j)&1) for i in range(p) for j in range(p)):continue
  if any(data['distances'][lift_coordinate(a,xs)]<model.d-min(a.bit_count(),model.m-a.bit_count()) for a in range(1,1<<p)):continue
  out.append(xs)
 return out

def physical(H,model,N):
 data=quotient_input(model,N);H,order,contexts=kernel_contexts(H,model.m,data['Gamma'],with_group=True)
 frame=standard(model.m);all_native=set();ours=set();records=[]
 for ci,c in enumerate(contexts):
  st,L=native(model,data,c);assert st['status']=='COMPLETE' and set(L)==set(oracle(model,data,c))
  all_native.update(L)
  for xs in L:
   a=solve(model,H,[lift_coordinate(x,data['units']) for x in xs],find_all=True)
   assert a['status']=='COMPLETE'
   for s in a['solutions']:
    S,D=model.construct(s['model_columns']);assert min(weight(w,model.n) for w in span(D) if w)>=model.d
    assert len(reduced_basis(s['standard_columns']))==len(data['Gamma'])
    ours.add(tuple(s['standard_columns']))
  path=HERE/'generic_group.txt';path.write_text(f'{len(c["group"])} {frame.p}\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
  gst,GL=native(model,data,c,group=path);assert gst['status']=='COMPLETE'
  canon={min(tuple(lift_coordinate(u,xs) for u in a) for a in c['group']) for xs in L}
  assert set(GL)==canon
  if L:assert native(model,data,c,fixed=L[0])[1]==[L[0]]
  records.append(dict(g=c['g'],t=c['t'],Z=len(L),Z_orbits=len(GL),prefix_group=len(c['group'])))
 # Independent full fixed-Z solver: enumerate all cost-valid independent Zs,
 # not just the leaves retained by the new kernel constraints.
 theirs=set();tested=0;p=frame.p;reps=[v for v,d in enumerate(data['distances']) if d>=model.d-1]
 for xs in itertools.product(reps,repeat=p):
  if len(reduced_basis(xs))!=p:continue
  if not z_condition(xs,data['Gamma'],len(data['Gamma'])):continue
  if any(data['distances'][lift_coordinate(a,xs)]<model.d-min(a.bit_count(),model.m-a.bit_count()) for a in range(1,1<<p)):continue
  full=[lift_coordinate(x,data['units']) for x in xs];answer=reference(model,frame,full,find_all=True);assert answer['status']=='COMPLETE';tested+=1
  for s in answer['solutions']:
   columns=s['adapted_columns'];mc=frame.model_columns(model,columns)
   if len(reduced_basis(columns))==len(data['Gamma']) and all(model.map_value(coordinates(h,frame.R+model.E)>>1,mc)==0 for h in H):
    assert xs in all_native,(H,xs,columns,contexts);theirs.add(tuple(columns))
 assert ours==theirs and ours
 return dict(parameters=[model.n,model.k,model.d],contexts=records,physical_maps=len(ours),reference_Z_cases=tested,site_order=order)

def main():
 start=time.monotonic();records=[]
 a=json.loads((HERE/'control.json').read_text());records.append(physical(a['H'],ExtensionModel(8,2,3,4,[240],a['C']),a['N']))
 S=json.loads((ROOT/'m8r1/fixed_radical_verification.json').read_text())['nonzero_radical_positive_control']['generators']
 H=[split(w,4,3)[0] for w in supported_subspace(dual(S,7),7,range(4))]
 records.append(physical(H,ExtensionModel(7,1,3,4,[240],[7]),[]))
 S=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'];D=dual(S,8)
 H=[split(w,4,4)[0] for w in supported_subspace(D,8,range(4))]
 C=[split(w,4,4)[1] for w in supported_subspace(D,8,range(4,8))]
 N=[split(w,4,4)[1] for w in supported_subspace(S,8,range(4,8))]
 records.append(physical(H,ExtensionModel(8,2,3,4,[240],C),N))
 # Artificial complete controls exercise both q3 cases, without a physical
 # distance claim. The pure-Z intersection and every matrix rank are checked.
 controls=[]
 class Toy:m=4;d=3
 for H,G in [([240,9,10,12],[0,0,0]),([240,9^32,10,12],[2,1,0])]:
  model=Toy();data=dict(labels=list(range(8)),distances=[0]+[3]*7,Gamma=G)
  H,order,cs=kernel_contexts(H,4,G,with_group=False)
  for c in cs:
   st,L=native(model,data,c);assert st['status']=='COMPLETE' and set(L)==set(oracle(model,data,c)) and L
   for limit in (0,1):assert native(model,data,c,limit=limit)[0]['status']=='LIMIT'
  controls.append(dict(g=cs[0]['g'],t=cs[0]['t'],identifications=len(cs),leaves_per_context=len(L)))
 report=dict(status='PASS',physical=records,artificial=controls,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'solve.py',HERE/'z_enum_fixed.c',HERE/'z_enum_fixed',ROOT/'model.py',ROOT/'mixed_frames/fixed_z.py',ROOT/'fixed_suffix_search.py']})
 atomic_json(HERE/'generic_verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
