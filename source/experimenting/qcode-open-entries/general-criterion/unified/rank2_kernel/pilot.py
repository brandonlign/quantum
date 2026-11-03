"""Bounded exact full-H target slices, independently isolated from stage4."""
import argparse,json,subprocess,sys,tempfile,time
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,span,weight
from rank2_kernel.common import local_R,contexts,quotient_input,payload,identification_orbits
from rank2_kernel.solve import solve
from rank2_kernel.symmetry import suffix_slice
from rank2_kernel.verify_z import BINARY
from model import lift_coordinate
from mixed_frames.stage3.common import atomic_json,sha

def search(model,H,frame,data,c,orbits,seconds):
 cmd=[str(BINARY),'--z-costs','--kernel-data','--gamma','--avoid-radical','--first-orbits']
 if seconds is not None:cmd+=['--seconds',str(seconds)]
 evidence=[];found=None
 with tempfile.TemporaryDirectory() as directory,tempfile.TemporaryFile('w+') as err:
  group=Path(directory)/'group.txt';group.write_text(f'{len(c["group"])} {frame.p}\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n');cmd+=['--prefix-group',str(group)]
  proc=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True)
  try:
   proc.stdin.write(payload(model,frame,data,c,orbits['first_representatives']));proc.stdin.close()
   for line in proc.stdout:
    xs=list(map(int,line.split()[1:]));answer=solve(model,H,frame,[lift_coordinate(x,data['units']) for x in xs],seconds=seconds);evidence.append(dict(X=xs,result=answer))
    if answer['status']=='FOUND':
     columns=answer['solutions'][0]['model_columns'];assert model.criterion(columns);S,D=model.construct(columns);minimum=min(weight(w,17) for w in span(D) if w);assert minimum>=6
     found=dict(S=S,D=D,minimum=minimum,columns=columns);proc.terminate();break
   code=proc.wait();err.seek(0)
   if found:stats=dict(status='FOUND_EARLY',leaves=len(evidence))
   else:assert code==0;stats=json.loads(err.read());assert stats['leaves']==len(evidence)
  finally:
   if proc.poll() is None:proc.terminate();proc.wait()
 complete=stats['status']=='COMPLETE' and all(e['result']['status']=='COMPLETE' and e['result']['fixed_z_excluded'] for e in evidence)
 return dict(status='FOUND' if found else 'EXCLUDED_IDENTIFICATION_SLICE' if complete else 'LIMIT',z=stats,leaves=evidence,found=found)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--classes',type=int,nargs='+',default=[0,12]);parser.add_argument('--max-N',type=int,default=1);parser.add_argument('--seconds',type=float,default=1);args=parser.parse_args()
 plan=json.loads((HERE/'suffix_plan.json').read_text());prefix=json.loads((HERE/'prefix.json').read_text())['prefix'];records=[]
 for row in plan['records']:
  if row['cls'] not in args.classes:continue
  model=ExtensionModel(17,3,6,8,local_R(8),row['C']);count=0
  for oi,o in enumerate(row['orbits']):
   if o['screen']:continue
   data=quotient_input(model,o['N'])
   for pi,p in enumerate(prefix[str(row['q'])]):
    H=p['H'];frame,cs=contexts(H,8,data['Gamma'],with_group=True);outer=identification_orbits(H,frame,cs)
    for orbit in outer['orbits']:
     ci=orbit['representative'];c=cs[ci];ss=suffix_slice(model,o['N'],frame,data,c);start=time.monotonic();r=search(model,H,frame,data,c,ss,args.seconds)
     r.update(cls=row['cls'],q=row['q'],N_orbit=oi,N=o['N'],H_index=pi,H=H,context_index=ci,context_orbit=orbit,context_orbits=len(outer['orbits']),prefix_group=len(c['group']),first_orbits=len(ss['first_representatives']),eligible=len(ss['eligible']),seconds=time.monotonic()-start)
     records.append(r);atomic_json(HERE/'pilot.json',dict(status='BOUNDED_PILOT',records=records,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'solve.py',HERE/'symmetry.py',BINARY,BINARY.with_suffix('.c'),HERE/'suffix_plan.json',HERE/'prefix.json']},scope='Bounded fixed-H,N,identification samples only; LIMIT is never an exclusion'))
     print(json.dumps({k:v for k,v in r.items() if k not in ('leaves','found','H','context_orbit')}),flush=True)
     if r['found']:return
   count+=1
   if count>=args.max_N:break
if __name__=='__main__':main()
