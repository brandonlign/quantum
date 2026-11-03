"""Bounded target benchmarks; a LIMIT is never an exclusion."""
import argparse,json,subprocess,sys,tempfile,time
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,gram,lift_coordinate,nullspace,reduced_basis,span,weight
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload,identification_orbits
from fixed_kernel.solve import solve
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha

def search(model,H,data,c,reps,seconds,*,group=True):
 cmd=[str(HERE/'z_enum_fixed'),'--gamma','--avoid-radical','--first-orbits','--kernel-data']
 if seconds is not None:cmd+=['--seconds',str(seconds)]
 tally=Counter();evidence=[];found=None
 with tempfile.TemporaryDirectory() as directory,tempfile.TemporaryFile('w+') as err:
  if group:
   p=Path(directory)/'prefix.txt';p.write_text(f'{len(c["group"])} {model.m-1}\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n');cmd+=['--prefix-group',str(p)]
  proc=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True)
  try:
   proc.stdin.write(payload(model,data,reps,c));proc.stdin.close()
   for line in proc.stdout:
    xs=list(map(int,line.split()[1:]));answer=solve(model,H,[lift_coordinate(x,data['units']) for x in xs],seconds=seconds)
    tally[answer['status']+':'+str(answer['reason'])]+=1
    evidence.append(dict(X=xs,result=answer))
    if answer['status']=='FOUND':found=answer;proc.terminate();break
   code=proc.wait();err.seek(0)
   if found:stats=dict(status='FOUND_EARLY',leaves=len(evidence))
   else:assert code==0;stats=json.loads(err.read());assert stats['leaves']==len(evidence)
  finally:
   if proc.poll() is None:proc.terminate();proc.wait()
 return dict(z=stats,tally=dict(tally),leaves=evidence,found=found)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=2);parser.add_argument('--count',type=int,default=1);args=parser.parse_args()
 source=ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json';rows=json.loads(source.read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records']
 out=[]
 for hull,q,t in [(1,1,0),(1,2,2),(3,2,0),(3,3,2),(5,3,0)]:
  chosen=[(i,r['representative']['generators']) for i,r in enumerate(rows) if r['representative']['hull']==hull][:args.count]
  for cls,C in chosen:
   model=ExtensionModel(17,3,6,8,[255<<8],C);rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),len(C))]
   N=[] if q==1 else (rad[:1] if q==2 else rad[:2]);data=quotient_input(model,N)
   assert len(reduced_basis(data['Gamma']))==14-2*q+t
   for ri,row in enumerate(profiles):
    if (row['q'],row['t'])!=(q,t):continue
    H=[normalize_word(h,row['p'],8) for h in row['H']];H,order,contexts=kernel_contexts(H,8,data['Gamma'],with_group=True);orbits=identification_orbits(H,8,contexts)
    ci=orbits['orbits'][0]['representative'];c=contexts[ci]
    eligible=[v for v,d in enumerate(data['distances']) if d>=5 and all(((v&r).bit_count()%2)==((c['tags'][0]>>j)&1) for j,r in enumerate(c['radical_coeff']))]
    start=time.monotonic();r=search(model,H,data,c,eligible,args.seconds)
    r.update(cls=cls,q=q,t=t,hull=hull,N=N,H=H,profile_index=ri,context=ci,identification_orbits=orbits,eligible_first=len(eligible),prefix_group=len(c['group']),seconds=time.monotonic()-start,scope='Bounded sample of ONE radical-identification slice; not a complete C,N or target exclusion')
    out.append(r);atomic_json(HERE/'pilot.json',dict(status='BOUNDED_PILOT',records=out,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'solve.py',HERE/'z_enum_fixed.c',HERE/'z_enum_fixed',source]}));print(json.dumps({k:v for k,v in r.items() if k not in ('leaves','found','H','identification_orbits')}),flush=True)
if __name__=='__main__':main()
