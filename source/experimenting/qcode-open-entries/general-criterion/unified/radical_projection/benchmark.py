"""Paired target-root pilots; no incomplete partition is an exclusion."""
import argparse,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,nullspace
from run_hull5_planes import quotient_input
from suffix_symmetries import syndrome_orbits
from rankone_light_budget.common import payload,OLD
from mixed_frames.stage3.common import atomic_json,sha
from radical_projection.verify import BINARY,valid

def measure(model,data,reps,mode,seconds):
 cmd=[str(OLD if mode=='old' else BINARY),'--gamma','--avoid-radical','--first-orbits','--seconds',str(seconds)]
 if mode=='projection':cmd+=['--projection']
 start=time.monotonic();p=subprocess.run(cmd,input=payload(model,data,reps),text=True,capture_output=True,check=True)
 xs=[list(map(int,l.split()[1:])) for l in p.stdout.splitlines()];stats=json.loads(p.stderr)
 assert stats['leaves']==len(xs)
 G=nullspace(data['Gamma'],len(data['Gamma']))
 passing=sum(valid(x,G,model.m,model.d) for x in xs)
 if mode=='projection':assert passing==len(xs)
 return dict(mode=mode,wall_seconds=time.monotonic()-start,stats=stats,passing_projection=passing,sample=xs[:10])

def main():
 p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=10);p.add_argument('--planes',type=int,default=1);args=p.parse_args()
 v=json.loads((HERE/'verification.json').read_text());assert v['status']=='PASS' and v['hashes'][str(BINARY.relative_to(ROOT))]==sha(BINARY)
 plan=json.loads((ROOT/'rankone_recenter/plane_plan.json').read_text());screen=json.loads((ROOT/'rankone_recenter/column_screen.json').read_text());lookup={r['cls']:r for r in screen['records']};results=[]
 for row in plan['records'][:3]:
  model=ExtensionModel(17,3,6,8,[255<<8],row['C'])
  indices=[i for i,o in enumerate(row['orbits']) if lookup[row['cls']]['orbits'][i]['reason'] is None][:args.planes]
  for i in indices:
   N=row['orbits'][i]['representative'];data=quotient_input(model,N);orbits=syndrome_orbits(model,data,N);reps=orbits['first_representatives'][:1];records=[]
   for mode in ('old','projection'):
    r=measure(model,data,reps,mode,args.seconds);records.append(r)
    print(json.dumps(dict(cls=row['cls'],orbit=i,first=reps,mode=mode,stats=r['stats'],passing_projection=r['passing_projection'])),flush=True)
   if all(r['stats']['status']=='COMPLETE' for r in records):assert records[0]['passing_projection']==records[1]['stats']['leaves']
   results.append(dict(cls=row['cls'],orbit=i,N=N,first=reps,measurements=records))
   atomic_json(HERE/'benchmark.json',dict(records=results,scope='Paired first-root pilots; LIMIT never excludes',backend_sha256=sha(BINARY)))
if __name__=='__main__':main()
