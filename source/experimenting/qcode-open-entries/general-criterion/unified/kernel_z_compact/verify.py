"""Byte-identical leaf streams AND node/prune counts against frozen native."""
import functools,json,random,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,dual,gram,lift_coordinate,nullspace,reduced_basis,split
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload
from gauge_fix import supported_subspace
from mixed_frames.stage3.common import atomic_json,sha
OLD=ROOT/'fixed_kernel/z_enum_fixed';NEW=HERE/'z_enum_compact'

def pair(model,data,c,*,group=False,fixed=(),limit=None,reps=None):
 if reps is None:reps=[v for v,d in enumerate(data['distances']) if d>=model.d-1]
 args=['--gamma','--avoid-radical','--first-orbits','--kernel-data']
 if fixed:args+=['--prefix',','.join(map(str,fixed))]
 if limit is not None:args+=['--max-nodes',str(limit)]
 with tempfile.TemporaryDirectory() as directory:
  if group:
   path=Path(directory)/'group.txt';path.write_text(f'{len(c["group"])} {model.m-1}\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n');args+=['--prefix-group',str(path)]
  out=[];text=payload(model,data,reps,c)
  for binary in [OLD,NEW]:
   start=time.monotonic();proc=subprocess.run([str(binary)]+args,input=text,capture_output=True,text=True,check=True,timeout=25);wall=time.monotonic()-start;stats=json.loads(proc.stderr)
   out.append((proc.stdout,stats,wall))
 assert out[0][0]==out[1][0]
 assert {k:v for k,v in out[0][1].items() if k!='seconds'}=={k:v for k,v in out[1][1].items() if k!='seconds'}
 return dict(status=out[0][1]['status'],nodes=out[0][1]['nodes'],leaves=out[0][1]['leaves'],old_seconds=out[0][1]['seconds'],new_seconds=out[1][1]['seconds'],old_wall=out[0][2],new_wall=out[1][2],N=len(data['labels']),columns=sum(d>=model.d-1 for d in data['distances']))

def main():
 start=time.monotonic();records=[];rng=random.Random(1713608)
 class Toy:m=4;d=3
 for g,t in [(0,2),(1,0),(2,0),(3,0),(1,2)]:
  b=6-g-t;G=[0]*b
  for i in range(0,b-g,2):G[i]=1<<(i+1);G[i+1]=1<<i
  for trial in range(8):
   data=dict(Gamma=G,labels=list(range(1<<b)),distances=[0]+[rng.randrange(1,5) for _ in range((1<<b)-1)])
   tags=[1<<j for j in range(g)]+[rng.randrange(1<<g) for _ in range(3-g)];tags.append(functools.reduce(int.__xor__,tags,0));F=[0]*3
   if t:F[g]=1<<(g+1);F[g+1]=1<<g
   c=dict(g=g,t=t,radical_coeff=nullspace(G,b),tags=tags,Gram=F)
   records.append(dict(kind='artificial',g=g,t=t,**pair(Toy(),data,c)))
 a=json.loads((ROOT/'fixed_kernel/control.json').read_text());fixtures=[(a['H'],ExtensionModel(8,2,3,4,[240],a['C']),a['N'])]
 S=json.loads((ROOT/'m8r1/fixed_radical_verification.json').read_text())['nonzero_radical_positive_control']['generators'];H=[split(w,4,3)[0] for w in supported_subspace(dual(S,7),7,range(4))];fixtures.append((H,ExtensionModel(7,1,3,4,[240],[7]),[]))
 S=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'];D=dual(S,8);H=[split(w,4,4)[0] for w in supported_subspace(D,8,range(4))];C=[split(w,4,4)[1] for w in supported_subspace(D,8,range(4,8))];N=[split(w,4,4)[1] for w in supported_subspace(S,8,range(4,8))];fixtures.append((H,ExtensionModel(8,2,3,4,[240],C),N))
 for H,model,N in fixtures:
  data=quotient_input(model,N);H,order,cs=kernel_contexts(H,4,data['Gamma'],with_group=True)
  for c in cs:
   for group in [False,True]:records.append(dict(kind='physical',parameters=[model.n,model.k,model.d],g=c['g'],group=group,**pair(model,data,c,group=group)))
 source=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
 from fixed_kernel.common import standard
 for control in json.loads((ROOT/'fixed_kernel/large_verification.json').read_text())['records']:
  C=source[control['cls']]['representative']['generators'];model=ExtensionModel(17,3,3,8,[255<<8],C);geometry=ExtensionModel(17,3,6,8,[255<<8],C);data=quotient_input(geometry,control['N']);H,order,cs=kernel_contexts(control['H'],8,data['Gamma']);frame=standard(8);full=frame.adapted_columns(model,control['known_columns']);small=[__import__('model').coordinates(x,data['units']) for x in full[:7]]
  for c in cs:
   if all((r&small[i]).bit_count()%2==((c['tags'][i]>>j)&1) for i in range(7) for j,r in enumerate(c['radical_coeff'])):break
  for pin in [7,6,5]:
   records.append(dict(kind='large_positive',g=c['g'],pin=pin,**pair(model,data,c,fixed=small[:pin])))
 report=dict(status='PASS',comparisons=len(records),records=records,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'PROOF.md',HERE/'z_enum_compact.c',NEW,OLD,ROOT/'fixed_kernel/z_enum_fixed.c',ROOT/'fixed_kernel/common.py',ROOT/'model.py']})
 atomic_json(HERE/'verification.json',report);print(json.dumps(dict(status='PASS',comparisons=len(records),seconds=report['seconds'])),flush=True)
if __name__=='__main__':main()
