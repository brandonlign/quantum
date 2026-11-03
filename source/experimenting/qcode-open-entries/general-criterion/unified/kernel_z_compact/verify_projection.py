"""Exact compact-column checks in the frozen radical-projection mode."""
import json,os,random,signal,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from radical_projection.verify import run as native
from model import ExtensionModel,dual,split
from gauge_fix import supported_subspace
from fixed_suffix_search import quotient_input
from m8r1.suffix_symmetries import syndrome_orbits
from rankone_light_budget.common import payload
from mixed_frames.stage3.common import atomic_json,sha
OLD=ROOT/'radical_projection/z_enum_projection';NEW=HERE/'z_enum_compact'

def main():
 start=time.monotonic();records=[];rng=random.Random(17366888)
 def compare(name,m,d,data):
  a,L=native(m,d,data,True,binary=OLD);b,K=native(m,d,data,True,binary=NEW)
  assert a['status']==b['status']=='COMPLETE' and L==K
  assert {k:v for k,v in a.items() if k!='seconds'}=={k:v for k,v in b.items() if k not in ('seconds','prefix_pruned')}
  records.append(dict(name=name,nodes=b['nodes'],leaves=b['leaves'],old_seconds=a['seconds'],new_seconds=b['seconds']))
 for b in range(3,7):
  for trial in range(4):
   pairs=(b-trial%4)//2;g=b-2*pairs
   if g>3:continue
   perm=list(range(b));rng.shuffle(perm);G=[0]*b
   for j in range(pairs):u,v=perm[2*j:2*j+2];G[u]=1<<v;G[v]=1<<u
   compare(f'artificial_b{b}_g{g}',4,3,dict(Gamma=G,distances=[0]+[rng.randrange(1,5) for _ in range((1<<b)-1)]))
 source=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'];C=[split(w,4,4)[1] for w in supported_subspace(dual(source,8),8,range(4,8))];N=[split(w,4,4)[1] for w in supported_subspace(source,8,range(4,8))]
 for model,N in [(ExtensionModel(5,0,3,3,[56],[]),[]),(ExtensionModel(7,1,3,4,[240],[7]),[]),(ExtensionModel(8,2,3,4,[240],C),N)]:compare(f'physical_{model.n}',model.m,model.d,quotient_input(model,N))
 pool=int((ROOT/'radical_projection/pool.pid').read_text());assert os.getpgid(pool)==pool
 assert '-m radical_projection.pool_run' in subprocess.check_output(['ps','-p',str(pool),'-o','command='],text=True)
 os.killpg(pool,signal.SIGSTOP)
 try:
  plan=json.loads((ROOT/'rankone_recenter/plane_plan.json').read_text())
  for cls,oi in [(150,0),(22354,96)]:
   row=next(r for r in plan['records'] if r['cls']==cls);N=row['orbits'][oi]['representative'];model=ExtensionModel(17,3,6,8,[65280],row['C']);data=quotient_input(model,N);reps=syndrome_orbits(model,data,N)['first_representatives'];text=payload(model,data,reps);out=[]
   for binary in [OLD,NEW]:
    t=time.monotonic();p=subprocess.run([str(binary),'--gamma','--avoid-radical','--first-orbits','--projection'],input=text,text=True,capture_output=True,check=True,timeout=50);out.append((p.stdout,json.loads(p.stderr),time.monotonic()-t))
   a,b=out;assert a[0]==b[0] and a[1]['status']==b[1]['status']=='COMPLETE'
   assert {k:v for k,v in a[1].items() if k!='seconds'}=={k:v for k,v in b[1].items() if k not in ('seconds','prefix_pruned')}
   records.append(dict(name=f'target_{cls}_{oi}',nodes=b[1]['nodes'],leaves=b[1]['leaves'],old_seconds=a[1]['seconds'],new_seconds=b[1]['seconds'],old_wall=a[2],new_wall=b[2],first_reps=len(reps),labels=len(data['labels']),columns=sum(d>=5 for d in data['distances'])))
 finally:os.killpg(pool,signal.SIGCONT)
 report=dict(status='PASS',records=records,seconds=time.monotonic()-start,scope='Complete exact leaf-stream and node-counter equality in projection mode; sampled target N-orbits only',hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'z_enum_compact.c',NEW,OLD,ROOT/'radical_projection/verify.py',ROOT/'radical_projection/verification.json',ROOT/'rankone_light_budget/common.py',ROOT/'fixed_suffix_search.py',ROOT/'model.py']})
 atomic_json(HERE/'projection_verification.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('records','hashes')}),flush=True)

if __name__=='__main__':main()
