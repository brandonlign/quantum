"""Exact recentered hull-five N-orbit sweep. LIMIT/ERROR never exclude."""
import argparse,json,os,sys,tempfile,time,subprocess
from collections import Counter
from functools import lru_cache
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,lift_coordinate,span,weight
from run_hull5_planes import quotient_input
from fixed_suffix_search import quotient_input as independent_quotient
from suffix_symmetries import syndrome_orbits
from signed_extension import solve_extension
from rankone_light_budget.common import payload
from mixed_frames.stage3.common import atomic_json,sha
from radical_projection.verify import BINARY,valid

@lru_cache(maxsize=8)
def target(C):return ExtensionModel(17,3,6,8,[255<<8],list(C))

def fingerprints():
 paths=[Path(__file__),BINARY,HERE/'z_enum_projection.c',HERE/'PROOF.md',HERE/'verification.json',ROOT/'rankone_recenter/plane_plan.json',ROOT/'rankone_recenter/column_screen.json',ROOT/'rankone_recenter/prepare.py',ROOT/'rankone_recenter/screen.py',ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json',ROOT/'m8r1/run_hull5_planes.py',ROOT/'m8r1/suffix_symmetries.py',ROOT/'fixed_suffix_search.py',ROOT/'model.py',ROOT/'signed_extension.py',ROOT/'two_site.py',ROOT/'kernel_bound.py']
 return {str(p.relative_to(ROOT)):sha(p) for p in paths}

def search(model,N,data,reps,seconds):
 cmd=[str(BINARY),'--gamma','--avoid-radical','--first-orbits','--projection']
 if seconds is not None:cmd+=['--seconds',str(seconds)]
 tally=Counter();evidence=[];found=None
 with tempfile.TemporaryFile('w+') as err:
  p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True)
  try:
   p.stdin.write(payload(model,data,reps));p.stdin.close()
   for line in p.stdout:
    fields=line.split();assert fields[0]=='LEAF' and len(fields)==model.m
    small=list(map(int,fields[1:]));xs=[lift_coordinate(x,data['units']) for x in small]
    answer=solve_extension(model,xs,minimum_stabilizer_weight=6)
    assert answer['status'] in ('FOUND','EXCLUDED','LIMIT')
    tally[answer['status']+':'+answer.get('reason','')]+=1
    evidence.append(dict(X=small,result=answer))
    if answer['status']=='FOUND':
     assert model.criterion(answer['model_columns']);S,D=model.construct(answer['model_columns'])
     assert min(weight(w,17) for w in span(D) if w)>=6
     found=dict(S=S,D=D,answer=answer);break
   if found:
    p.terminate();p.wait();stats=dict(status='FOUND_EARLY',leaves=len(evidence))
   else:
    assert p.wait()==0;err.seek(0);stats=json.loads(err.read())
    assert stats['status'] in ('COMPLETE','LIMIT') and stats['leaves']==len(evidence)
  finally:
   if p.poll() is None:p.terminate();p.wait()
 exhausted=stats['status']=='COMPLETE' and all(e['result']['status']=='EXCLUDED' for e in evidence)
 return dict(status='FOUND' if found else ('EXCLUDED_FIXED_N' if exhausted else 'LIMIT'),z=stats,tally=dict(tally),leaves=evidence,found=found)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float);parser.add_argument('--max-tasks',type=int);args=parser.parse_args();assert __debug__
 verified=json.loads((HERE/'verification.json').read_text());assert verified['status']=='PASS'
 for path,value in verified['hashes'].items():assert sha(ROOT/path)==value,path
 hashes=fingerprints();plan=json.loads((ROOT/'rankone_recenter/plane_plan.json').read_text());screen=json.loads((ROOT/'rankone_recenter/column_screen.json').read_text());lookup={r['cls']:r for r in screen['records']};count=0
 assert plan['status']=='PASS' and screen['status']=='PASS'
 assert plan['source_sha256']==sha(ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json')
 assert plan['preparation_sha256']==sha(ROOT/'rankone_recenter/prepare.py')
 (HERE/'runs').mkdir(exist_ok=True)
 for row in plan['records']:
  model=target(tuple(row['C']));assert model.hull==5 and model.C_distance>=6
  for i,orbit in enumerate(row['orbits']):
   if lookup[row['cls']]['orbits'][i]['reason'] is not None:continue
   path=HERE/'runs'/f'{row["cls"]}_{i}.json'
   if path.exists():
    previous=json.loads(path.read_text());assert previous['hashes']==hashes
    if previous['status']=='EXCLUDED_FIXED_N':continue
    if previous['status']=='FOUND':raise RuntimeError('Inspect saved pure-code candidate')
   start=time.monotonic();N=orbit['representative'];data=quotient_input(model,N)
   check=independent_quotient(model,N)
   assert all(data[key]==check[key] for key in ('units','labels','distances','Gamma'))
   orbits=syndrome_orbits(model,data,N)
   result=search(model,N,data,orbits['first_representatives'],args.seconds)
   result.update(cls=row['cls'],orbit_index=i,N=N,orbit_data=orbits,hashes=hashes,seconds=time.monotonic()-start,scope='One exact-N orbit in hull-five exclusive rank-one column, conditional on suffix classification')
   atomic_json(path,result);print(json.dumps({k:result[k] for k in ('cls','orbit_index','status','seconds','z','tally')}),flush=True);count+=1
   if result['status']=='FOUND' or (args.max_tasks is not None and count>=args.max_tasks):return
 atomic_json(HERE/'sweep_finished.json',dict(status='FINISHED_REQUIRES_AGGREGATION',hashes=hashes))
if __name__=='__main__':main()
