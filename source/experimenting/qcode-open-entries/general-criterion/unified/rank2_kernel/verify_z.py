"""Independent generic-cost oracle, complete small physical maps and positives."""
import itertools,json,random,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,lift_coordinate,nullspace,reduced_basis
from rank2_kernel.common import local_R,contexts,quotient_input,payload
from rank2_kernel.solve import solve
from mixed_frames.fixed_z import solve as reference
from mixed_frames.stage3.common import atomic_json,sha
BINARY=ROOT/'rank2_kernel_backend/z_enum_cost'

def run(model,frame,data,c,*,fixed=(),group=False,limit=None):
 reps=[v for v,d in enumerate(data['distances']) if d>=model.d-frame.costs[1]]
 cmd=[str(BINARY),'--z-costs','--kernel-data','--gamma','--avoid-radical','--first-orbits']
 if fixed:cmd+=['--prefix',','.join(map(str,fixed))]
 if limit is not None:cmd+=['--max-nodes',str(limit)]
 if group:
  path=HERE/'test_group.txt';path.write_text(f'{len(c["group"])} {frame.p}\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n');cmd+=['--prefix-group',str(path)]
 proc=subprocess.run(cmd,input=payload(model,frame,data,c,reps),text=True,capture_output=True,timeout=20)
 assert proc.returncode==0,(proc.stderr,frame.p,data['Gamma'],c)
 stats=json.loads(proc.stderr);L=[tuple(map(int,l.split()[1:])) for l in proc.stdout.splitlines()];assert len(L)==len(set(L))==stats['leaves']
 return stats,L

def oracle(model,frame,data,c,fixed=()):
 p=frame.p;g=c['g'];pre={lift_coordinate(a,data['Gamma']):a for a in range(len(data['labels']))};domains=[]
 for i in range(p):
  domains.append([v for v,d in enumerate(data['distances']) if (i>=len(fixed) or v==fixed[i]) and d>=model.d-frame.costs[1<<i] and all((r&v).bit_count()%2==((c['tags'][i]>>j)&1) for j,r in enumerate(c['radical_coeff']))])
 out=[]
 for xs in itertools.product(*domains):
  if len(reduced_basis(xs))!=p:continue
  ts=[xs[i]^lift_coordinate(coordinates(c['tags'][i],c['tags'][:g]),xs[:g]) for i in range(p)]
  assert all(t in pre for t in ts)
  if any((pre[ts[i]]&ts[j]).bit_count()%2!=((c['Gram'][i]>>j)&1) for i in range(p) for j in range(p)):continue
  if any(frame.costs[a]+data['distances'][lift_coordinate(a,xs)]<model.d for a in range(1,1<<p)):continue
  out.append(xs)
 return out

def main():
 start=time.monotonic();rng=random.Random(88817136);art=[]
 class Object:pass
 for p,g,t in [(2,0,0),(2,1,0),(2,2,0),(3,3,0),(3,0,2),(3,1,2)]:
  b=2*p-g-t;G=[0]*b
  for i in range(0,b-g,2):G[i]=1<<(i+1);G[i+1]=1<<i
  for trial in range(4):
   model=Object();model.d=3;frame=Object();frame.p=p;frame.costs=[0]+([2]*((1<<p)-1) if trial==0 else [rng.randrange(1,5) for _ in range((1<<p)-1)])
   data=dict(Gamma=G,labels=list(range(1<<b)),distances=[0]+([1]*((1<<b)-1) if trial==0 else [rng.randrange(1,4) for _ in range((1<<b)-1)]))
   tags=[1<<i for i in range(g)]+[rng.randrange(1<<g) for _ in range(p-g)];tags.append(__import__('functools').reduce(int.__xor__,tags,0));F=[0]*p
   if t:F[g]=1<<(g+1);F[g+1]=1<<g
   c=dict(g=g,t=t,radical_coeff=nullspace(G,b),tags=tags,Gram=F)
   st,L=run(model,frame,data,c);assert st['status']=='COMPLETE' and set(L)==set(oracle(model,frame,data,c))
   assert run(model,frame,data,c,limit=0)[0]['status']=='LIMIT'
   if L:assert run(model,frame,data,c,limit=1)[0]['status']=='LIMIT'
   art.append(dict(p=p,g=g,t=t,trial=trial,Z=len(L)))
 controls=[]
 for row in json.loads((HERE/'controls.json').read_text())['records']:
  n,k,d=row['parameters'];model=ExtensionModel(n,k,d,row['m'],local_R(row['m']),row['C']);data=quotient_input(model,row['N']);frame,cs=contexts(row['H'],row['m'],data['Gamma'],with_group=n<=9);c=cs[row['context_index']];known=tuple(row['small_columns'][:frame.p]);fixed=() if n<=9 else known[:-2]
  st,L=run(model,frame,data,c,fixed=fixed);assert st['status']=='COMPLETE' and known in L and set(L)==set(oracle(model,frame,data,c,fixed))
  group_count=None;maps=0
  if n<=9:
   gst,GL=run(model,frame,data,c,group=True);assert gst['status']=='COMPLETE'
   assert set(GL)=={min(tuple(lift_coordinate(a,xs) for a in U) for U in c['group']) for xs in L};group_count=len(GL)
   for xs in L:
    X=[lift_coordinate(x,data['units']) for x in xs];answer=solve(model,row['H'],frame,X,find_all=True);ref=reference(model,frame,X,find_all=True)
    assert answer['status']==ref['status']=='COMPLETE'
    mine={tuple(s['model_columns']) for s in answer['solutions']};theirs={tuple(s['model_columns']) for s in ref['solutions'] if len(reduced_basis(s['model_columns']))==len(data['Gamma']) and all(model.map_value(coordinates(h,frame.R+model.E)>>2,s['model_columns'])==0 for h in row['H'])};assert mine==theirs;maps+=len(mine)
  X=[lift_coordinate(x,data['units']) for x in known];answer=solve(model,row['H'],frame,X);assert answer['status']=='FOUND'
  controls.append(dict(parameters=row['parameters'],q=c['g'],Z=len(L),Z_orbits=group_count,complete_physical_maps=maps,known_positive_retained=True,pinned=len(fixed)))
 report=dict(status='PASS',artificial=art,physical=controls,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),BINARY,BINARY.with_suffix('.c'),HERE/'common.py',HERE/'solve.py',HERE/'controls.json',ROOT/'model.py',ROOT/'mixed_frames/fixed_z.py']})
 atomic_json(HERE/'z_verification.json',report);print(json.dumps(dict(status='PASS',artificial_cases=len(art),physical=controls,seconds=report['seconds'])),flush=True)
if __name__=='__main__':main()
